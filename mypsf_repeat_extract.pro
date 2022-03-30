PRO CALEVAL, X, P, YMOD

  YMOD = P[0] + P[1]*X
  YMOD = image_shift(YMOD,P[2],P[3])

END

;---------------------

; Contrary to StarFinder PSF_REPEAT_EXTRACT
; this routine takes into account that the references
; stars can also be secondary sources that must be subtracted.
; This is helpful in very dense fields when the reference stars are
; closely spaced.
;
; KEYWORDS
;
; MINDIST    Secondary stars at less than MINDIST separation of a PSF
;            reference star will not be subtracted. This is an
;            untested feature, but I suspect that it avoids rarely
;            occurring problems with PSF extraction in very crowded
;            fields.
;            The default value of MINDIST is 3 pixels.
; LAST MODIFICATION
;
; Introduced MINDIST keyword to avoid that stars very close to the
; reference star (or the reference star itself) can be subtracted as secondaries.
; 
; Activated the CALEVAL procedure again that scales the photometry of
; secondaries to the current best-estimate of the PSF.
; Rainer Schoedel 19 Oct 2021

PRO MYPSF_REPEAT_EXTRACT, x_psf, y_psf, x_stars, y_stars, f_stars, image, psf, normrad, REBFAC=rebfac, N_SIGMA = n_sigma, REL_THRESH = rel_thresh, BACKGROUND=background, CORE = core, DEBUG = debug, ITER = iter, MINDIST = mindist, PSF_NOISE = psf_noise, MASKRAD = mrad, UNWEIGHTED=unweighted

 if not(keyword_set(debug)) then debug = 0
 if not(keyword_set(rebfac)) then rebfac = 1
 if not(keyword_set(mindist)) then mindist = 2.


 sz = size(psf)
 boxsize = sz[1]
 boxhw = boxsize/2

 sz = size(image)
 nax1 = sz[1]
 nax2= sz[2]

 ; apply rebin factor
 ; use new variabels for rebinned quantities so that the parameters
 ; in the callin gprogram will not be changed
 nax1 = rebfac*nax1
 nax2 = rebfac*nax2
 normrad = rebfac * normrad
 mindist = rebfac * mindist
 boxhw = rebfac * boxhw
 boxsize = rebfac*boxsize
 x_ref = rebfac * x_psf
 y_ref = rebfac * y_psf
 x_s = rebfac * x_stars
 y_s = rebfac * y_stars
 f_s = f_stars

 nref = n_elements(x_ref)
 xint = round(x_ref)
 yint = round(y_ref)

psfim = image
if (rebfac gt 1) then begin
 psfim = CREBIN(psfim,nax1,nax2,/TOTAL)
 psf = CREBIN(psf,boxsize,boxsize,/TOTAL)
endif
psf = psf/total(psf)

if KEYWORD_SET(background) then begin
 bg = CREBIN(background,nax1,nax2,/TOTAL)
 psfim = psfim - bg
endif

bg_dummy = replicate(1,boxsize,boxsize)
bg_dummy = CIRC_MASK(bg_dummy,boxhw,boxhw,mrad,/INNER)

; loop over iterations
for it = 1, iter do begin

 psf = CIRC_MASK(psf,boxhw,boxhw,mrad,BORDER=0)
 psf = psf/total(psf)

 ; extract sub-images centered on reference stars
 ; and subtract secondary sources
 ;-------------------------------------------------
 stack = fltarr(boxsize,boxsize,nref)
 psfmasks = fltarr(boxsize,boxsize,nref)
 for iref = 0, nref-1 do begin
  print, xint[iref], yint[iref]
  sub_arrays, psfim, xint[iref], yint[iref], boxsize, slice, masks
  psfmasks[*,*,iref] = masks
  distances = sqrt((x_s-x_ref[iref])^2 + (y_s-y_ref[iref])^2)
  in_box = where(distances lt boxsize,in_count)
  other = where(distances lt boxsize and distances gt mindist)
  if (other[0] gt -1) then begin
   x_in_box = x_s[in_box] - (xint[iref] - boxhw)
   y_in_box = y_s[in_box] - (yint[iref] - boxhw)
   f_in_box = f_s[in_box]
   model = image_model(x_in_box,y_in_box,f_in_box,boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
   P = [0.,1.,0.,0.]
   W = model
   res = mpcurvefit(model,slice,W,P,sigma,FUNCTION_NAME='CALEVAL',/NODERIVATIVE,/QUIET)
   x_other = x_s[other] - (xint[iref] - boxhw)
   y_other = y_s[other] - (yint[iref] - boxhw)
   f_other = f_s[other]
   other_stars = image_model(x_other,y_other,f_other,boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
 ;  if (debug gt 0) then begin
 ;   print, P
 ;   writefits, 'rawslice.fits', slice
 ;   writefits, 'other_stars.fits', other_stars
 ;   writefits, 'diff.fits', slice - other_stars
 ;   writefits, 'other_stars.fits', image_shift(other_stars*P[1],P[2],P[3])
 ;   writefits, 'diff_fit.fits', slice - image_shift(other_stars*P[1],P[2],P[3])
 ;   STOP
 ;  endif
   slice = slice - image_shift(P[0] + other_stars*P[1],P[2],P[3])
  endif
  stack[*,*,iref] = slice
 endfor

; sub-pixel shift, background subtract (is essential!), and normalize the slices of the stack
psf_weights = intarr(nref) ; for  weighting
for iref = 0, nref-1 do begin
 subim = stack[*,*,iref]
 mask = psfmasks[*,*,iref]
 bgind = where(bg_dummy gt 0 and mask gt 0)
 subim = subim - median(subim[bgind])
; Sub-pixel centering
  subim = centroider(subim, XC = boxhw, YC = boxhw, _EXTRA = extra)
  normim = subim
  if min(normim lt 0) then normim[where(normim lt 0)] = 0
  normim = circ_mask(normim, boxhw, boxhw, normrad)
  ; normalization will not be correct for saturated stars
  fnorm = total(normim)
  psf_weights[iref] = round(sqrt(fnorm))
  subim = subim/fnorm
  stack[*,*,iref] = subim
endfor

; Define weights for median superposition
;if  not keyword_set(unweighted)  then begin
;   w = fltarr(nref)
;   for  n = 0L, nref - 1  do  w[n] = sqrt(stack[boxhw,boxhw,n] > 0)
;endif 

if (debug gt 0) then begin
 writefits, 'stack.fits', stack
endif

; create PSF through stack median
;psf = stack_deconv(stack, MASK=psfmasks)
;STOP
if keyword_set(unweighted) then psf_weights[*] = 1
psf = stack_median(stack, WEIGHTS=psf_weights, MASK=psfmasks)
psf_noise = stack_error(stack, WEIGHTS=psf_weights, MASK=psfmasks)
;writefits, 'hhh.fits', psf
;STOP

 ; OPTIONAL: threshold cut PSF and extract only contiguous signal
 ; ------------------------------------------------------------------

 IF KEYWORD_SET(n_sigma) and KEYWORD_SET(mrad) then begin
   if (debug gt 0) then begin
      ;writefits, 'bgring.fits', psf*bgring
      writefits, 'psf_raw.fits', psf
   endif
   if (rel_thresh eq 1) then begin
    ; ringmask for background/noise estimation
    dummy = fltarr(boxsize,boxsize)
    dummy[*,*] = 1
    bgring = circ_mask(dummy,boxhw,boxhw,mrad,BORDER=0, INNER=1)
    bgind = where(bgring gt 0)
    print, 'PSF background: ' + string(median(psf[bgind]))
    print, 'PSF std dev: ' + string(stddev(psf[bgind]))
    if nref ge 5 then begin
     noise = psf_noise
     threshold =  median(psf[bgind]) + n_sigma * noise 
     threshold = median(threshold) ; necessary because threshold is an array in this case
    endif else begin
     noise = stddev(psf[bgind])
     threshold =  median(psf[bgind]) + n_sigma * noise
    endelse 
  endif else begin
   threshold = n_sigma * max(psf)
  endelse
  if (CORE gt 0 ) then begin
   psf = image_core(psf,threshold,/SUBTRACT)
  endif else begin
   psf = psf - threshold
   suppress = where(psf lt 0, count)        
   if (count gt 0) then psf[suppress] = 0
  endelse
 endif

 ; mask PSF for all but last iteration
 if (it lt iter) then begin
  MASK_PSF, psf, mrad, PSF_MASKED=psf_masked, PSF_OFFSET=psf_offset, WINGS=wings
  psf = psf_masked
  neg = where(psf lt 0)
  if (neg[0] gt -1) then psf[neg] = 0
 endif
 psf_norm = total(psf)
 psf_noise = psf_noise/psf_norm
 psf = psf/psf_norm          ; normalization of PSF
 
; Final centroiding: not necessary?
; psfcen = centroid(circ_mask(psf,boxhw,boxhw,mrad,BORDER=0))   
 if (debug gt 0) then begin
  writefits, 'psf_it'+strn(it)+'.fits', psf
 endif

endfor ; end loop over iter

if (rebfac gt 1) then begin
 boxsize = boxsize/rebfac
 psf = CREBIN(psf,boxsize,boxsize,/TOTAL)
endif
psfnorm = total(psf)
psf = psf/psfnorm
psf_noise = psf_noise/psfnorm

END
