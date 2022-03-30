;PRO CALEVAL, X, P, YMOD

;  YMOD = P[0] + P[1]*X
;  YMOD = image_shift(YMOD,P[2],P[3])

;END


;---------------------

; Contrary to StarFinder PSF_REPEAT_EXTRACT
; this routine takes into account that the references
; stars can also be secondary sources that must be subtracted.
; This is helpful in very dense fields when the reference stars are
; closely spaced.


PRO PSF_EXTRACT_REPEAT, x_psf, y_psf, x_stars, y_stars, f_stars, image, psf, NORM_MAX = norm_max, REBFAC=rebfac, N_SIGMA = n_sigma, REL_THRESH = rel_thresh, MASKRAD=maskrad, BACKGROUND=background, USE_BG = use_bg, CORE = core, DEBUG = debug, ITER = iter, CLEAN_STACK = clean_stack, PSFMASKS = psfmasks


 dmax = 2. * rebfac

 if not(keyword_set(debug)) then debug = 0
 if not(keyword_set(norm_max)) then norm_max = 0

 sz = size(psf)
 boxsize = sz[1]
 boxhw = boxsize/2 

 sz = size(image)
 nax1 = sz[1]
 nax2= sz[2]

 ; apply rebin factor
 ; use new variabels for rebinned quantities so that the parameters
 ; in the calling program will not be changed
 if keyword_set(maskrad) then mrad = rebfac*maskrad
 nax1 = rebfac*nax1
 nax2 = rebfac*nax2
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

 psfim = FREBIN(image,nax1,nax2,/TOTAL)
 psf = FREBIN(psf,boxsize,boxsize,/TOTAL)
 psf = psf/total(psf)

if KEYWORD_SET(background) and (use_bg gt 0) then begin
 bg = FREBIN(background,nax1,nax2,/TOTAL)
 psfim = psfim - bg
endif

; loop over iterations
for it = 1, iter do begin

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
 nearby = where(distances lt boxsize)
 if (nearby[0] gt -1) then begin
  x_near = x_s[nearby] - (xint[iref] - boxhw)
  y_near = y_s[nearby] - (yint[iref] - boxhw)
  f_near = f_s[nearby]
;  model = image_model(x_near,y_near,f_near,boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
;  P = [0.,1.,0.,0.]
;  W = model
;  res = mpcurvefit(model,slice,W,P,sigma,FUNCTION_NAME='CALEVAL',/NODERIVATIVE,/QUIET)
  compare_lists, x_near, y_near, x_ref[iref]-(xint[iref]-boxhw), y_ref[iref]-(yint[iref]-boxhw), x1c, y1c, x2c ,y2c, SUB1=other_stars, MAX_DISTANCE=dmax
  if (other_stars[0] gt -1) then begin
   other_stars = image_model(x_near[other_stars],y_near[other_stars],f_near[other_stars],boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
 ;  if (debug gt 0) then begin
 ;   print, P
 ;   writefits, 'rawslice.fits', slice
 ;   writefits, 'other_stars.fits', image_shift(other_stars*P[1],P[2],P[3])
 ;   writefits, 'diff.fits', slice - image_shift(other_stars*P[1],P[2],P[3])
 ;   STOP
 ;  endif
   slice = slice - other_stars ; image_shift(other_stars*P[1],P[2],P[3])
  endif
 endif
 stack[*,*,iref] = slice
endfor

; sub-pixel shift and normalize the slices of the stack
for iref = 0, nref-1 do begin
 subim = stack[*,*,iref]
; Sub-pixel centering
  subim = centroider(subim, XC = boxhw, YC = boxhw, _EXTRA = extra)
 if (norm_max eq 1) then begin
  subim = subim/max(subim)
 endif else begin
  nrad = round(2*fwhm(psf))
  normim = circ_mask(subim, boxhw, boxhw, nrad,BORD=0)
  subim = subim/total(normim)
 endelse
 stack[*,*,iref] = subim
endfor
if (debug gt 0) then begin
 writefits, 'stack.fits', stack
endif
clean_stack = stack

; create PSF through stack median 
psf = stack_median(stack, MASK=psfmasks)


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
    noise = stddev(psf[bgind])
    threshold =  median(psf[bgind]) + n_sigma * noise
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
; Final centroiding: not necessary?
; psfcen = centroid(circ_mask(psf,boxhw,boxhw,mrad,BORDER=0))   
 psf = psf/total(psf)
 if (debug gt 0) then begin
  writefits, 'psf_core.fits', psf
 endif

endfor ; end loop over iter

END
