PRO CALEVAL, X, P, YMOD

  YMOD = P[0] + P[1]*X
  YMOD = image_shift(YMOD,P[2],P[3])

END


;---------------------

PRO HOLO_NOSUB, indir, inlist, outnam,  maskrad, normrad, nsigma, rebfac, refsources, starlist, maskdir=maskdir, masklist=masklist, DEBUG = debug, iter=iter, AIRY=airy, BORDER=bord, OUT_ITER=out_iter, PSFOUT = psfout, UNWEIGHTED = unweighted, SUBPIX=subpix, PSFAVG=psfavg, CLIP=clip, tmpdir=tmpdir, BOXHW=boxhw, NSUB=nsub, MINSUPP=minsupp, MAXSUPP = maxsupp, WEIGHTFRAMES=weightframes

; VERSION DATE: 8 September 2012


withmask = KEYWORD_SET(masklist)
if (not KEYWORD_SET(minsupp)) then minsupp = 0
if (not KEYWORD_SET(maxsupp)) then maxsupp = 1.5

; Identical to holo7sub.pro except that it 
; subtracts the noise threshold from the final PSF
; -------------------------------------------------

 ; star list
  readcol, starlist, xx, yy, ff
  xx = rebfac * xx
  yy = rebfac * yy
  ff = rebfac^2 * ff
; read list of input cubes
name = ''
openr, inp, inlist, /get_lun
ncubes = 0L
while (not EOF(inp)) do begin
  readf, inp, name
  ncubes = ncubes + 1
endwhile
list = strarr(ncubes)
point_lun, inp, 0
for ic = 0, ncubes -1 do begin
  readf, inp, name
  list[ic] = name
endfor
free_lun, inp
print
print, 'input file list: ' 
print, list
print

; read list of  FOV masks
if (withmask) then begin
 mlist = strarr(ncubes)
 openr, inp, masklist, /get_lun
 for ic = 0, ncubes -1 do begin
  readf, inp, name
  mlist[ic] = name
 endfor
 free_lun, inp
 print
 print, 'mask list: ' 
 print, mlist
 print
endif

im = readfits(indir + list[0], NSLICE=0)
sz = size(im)
nax1 = sz[1]
nax2 = sz[2]

; multiply with rebin factor
maskrad = rebfac*maskrad
nax1 = rebfac*nax1
nax2 = rebfac*nax2
normrad = rebfac*normrad
boxhw = rebfac*boxhw
; pixel indixes for subarray extraction
boxsize = long(2*boxhw + 1)

cen1 = nax1/2
cen2 = nax2/2
center = [cen1,cen2]

; dummy needed  ringmask for background estimation
dummy = fltarr(boxsize,boxsize)
dummy[*,*] = 1
bgring = circ_mask(dummy,boxhw,boxhw,maskrad,BORDER=0, INNER=1)
bgind = where(bgring gt 0)

; apodization with telescope OTF
; be careful to provide an adequately sampled PSF 
; of the correct size 
; (re-binning of PSF can be done but is non-optimal)
airy = airy/total(airy)
APOD = FFT(airy)*(nax1*nax2)
OTF = ABS(APOD)

finumer = fltarr(nax1,nax2)
fidenom = fltarr(nax1,nax2)
numersubreal = fltarr(nax1,nax2,nsub)
denomsubreal = fltarr(nax1,nax2,nsub)
numersubim = fltarr(nax1,nax2,nsub)
denomsubim = fltarr(nax1,nax2,nsub)
ngood = 0L

if (withmask) then begin
 fimasknumer = fltarr(nax1,nax2)
 masknumersubreal = fltarr(nax1,nax2,nsub)
 masknumersubim = fltarr(nax1,nax2,nsub)
endif

; reference sources
x_psf = refsources[0,*]
y_psf = refsources[1,*]
x_psf = rebfac * x_psf
y_psf = rebfac * y_psf
nref = n_elements(x_psf)
xint = round(x_psf)
yint = round(y_psf)


for ic = 0, ncubes-1 do begin  ; start loop over all input cubes

  ; initialize variables for holography
  ;---------------------------------------

  masknumer = fltarr(nax1,nax2)
  numer = fltarr(nax1,nax2)
  numer[*,*] = 0.0
  denom = fltarr(nax1,nax2)
  denom[*,*] = 0.0
  nthis = 0L

  ; read speckle data cube
  ; --------------------------
  if (withmask) then mcube = readfits(maskdir + mlist[ic])
  cube = readfits(indir + list[ic], header)
  cube = cube
  sz = size(cube)
  if (sz[0] gt 2) then  nim = sz[3] else nim = 1
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  print, 'Read in cube: ' + list[ic]
  print, 'Number of frames in this cube: ' + strtrim(string(nim),2)
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'

  rcims = fltarr(nax1,nax2,nim)

  if keyword_set(psfout) then psfs = fltarr(nax1,nax2,nim)
  ngpsf = 0L

  for i = 0L, nim-1 do begin  ; start loop over frames in this cube

    ; load ith image 
    im = cube[*,*,i]

   if (withmask) then begin
    thisfov = mcube[*,*,i]  
;    zeroind = where(thisfov lt 1, complement = ones)
    ; fov can have pixels with values << 1
    ; set them to zero to avoid problems with routines that
    ; use the masks
;    if (zeroind[0] gt -1) then thisfov[zeroind] = 0
;    thisfov[ones] = 1
    im = im * thisfov ; in case im has non-zero values outside of FOV and for correct weighting of mask edges
   endif
    
   ; rebin image and FOV mask and adjust stellar positions
    if (rebfac gt 1.0) then begin
  ;  writefits, tmpdir + 'origim.fits', im
     im = FREBIN(im,nax1,nax2,/TOTAL)
     if (withmask) then begin
      thisfov = FREBIN(thisfov,nax1,nax2) ; Don normalize rebinned FOV!
     endif
    endif
 
  ; PSF extraction
  ; ################

      psfim = im
      sub_arrays, psfim, xint, yint, boxsize, stack, masks
      ; modify psfmasks if a reference star may not be contained in all
      ; cubes because of ditherring
      ; ---------------------------
      if withmask then begin
       sub_arrays, thisfov, xint, yint, boxsize, psfmasks, masks
       for iref = 0, nref-1 do masks[*,*,iref] = masks[*,*,iref]*psfmasks[*,*,iref]
      endif
      ; save non-normalized stack for photometry
      rawstack = stack
      w_frame = 0.0 ; initialize frame weight
      for iref = 0, nref-1 do begin
       subim = stack[*,*,iref]
       if (subpix gt 0) then begin
        subim = image_shift(subim,xint[iref]-x_psf[iref],yint[iref]-y_psf[iref])
       endif
       normim = subim - median(subim[bgind])
       ; the next IF statement avoids a problem
       ; if the entire subim = 0 (can occur when a reference
       ; star lies outside the crrent FOV
       if min(normim lt 0) then normim[where(normim lt 0)] = 0
       normim = circ_mask(normim, boxhw, boxhw, normrad,BORD=0)
       w_frame = w_frame + max(normim)
       if (unweighted eq 1 and total(normim) gt 0) then begin
        subim = subim/total(normim)
       endif
       stack[*,*,iref] = subim
      endfor
      w_frame = 1.e-3 * w_frame/nref
;      print, ' Frame weight: ' + strn(w_frame)
      psf = stack_avg(stack, MASK=masks, AVGYN=0)
 
       ; mask PSF
      ; -----------------------
      if (debug gt 0) then begin
       writefits, tmpdir + 'im.fits', im
       if withmask then writefits, tmpdir + 'thisfov.fits', thisfov
       writefits, tmpdir + 'bgring_0.fits', psf*bgring
       writefits, tmpdir + 'psf_0raw.fits', psf
       writefits, tmpdir + 'stack_0.fits', stack*masks
      endif
;      print, 'Background: ' + string(median(psf[bgind]))
      psf = psf - median(psf[bgind])  ; subtract background
      noise = stddev(psf[bgind])
      suppress = where(psf lt nsigma[0] * noise, count)        
      if (count gt 0) then psf[suppress] = 0
      psfcen = centroid(circ_mask(psf,boxhw,boxhw,maskrad,BORDER=0))   
      psf = circ_mask(psf,psfcen[0],psfcen[1],maskrad,BORDER=0)
      psf = psf/total(psf)
      if (debug gt 0) then begin
        writefits, tmpdir + 'psf_0.fits', psf
      endif

   
  ; improval of PSF with known sources
  ; ####################################

      for it = 0, iter-1 do begin

       ; clean surroundings of reference stars from secondary sources
       ; ------------------------------------------------------------
      stack = fltarr(boxsize,boxsize,nref)
      ; use psfmasks if a reference star may not be contained in all
      ; cubes because of ditherring
      ; ---------------------------
      psfmasks = fltarr(boxsize,boxsize,nref)
      for iref = 0, nref-1 do begin
       sub_arrays, im, xint[iref], yint[iref], boxsize, slice, masks
       if withmask then begin
        sub_arrays, thisfov, xint[iref], yint[iref], boxsize, maskslice, masks
        psfmasks[*,*,iref] = maskslice * masks
       endif  else psfmasks[*,*,iref] = masks
       distances = sqrt((xx-x_psf[iref])^2 + (yy-y_psf[iref])^2)
       nearby = where(distances lt boxsize)
       if (nearby[0] gt -1) then begin
        x_near = xx[nearby] - (xint[iref] - boxhw)
        y_near = yy[nearby] - (yint[iref] - boxhw)
        f_near = ff[nearby]
        model = image_model(x_near,y_near,f_near,boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
        P = [0.,1.,0.,0.]
        W = model
        res = mpcurvefit(model,slice,W,P,sigma,FUNCTION_NAME='CALEVAL',/NODERIVATIVE,/QUIET)
        compare_lists, x_near, y_near, x_psf[iref]-(xint[iref]-boxhw), y_psf[iref]-(yint[iref]-boxhw), x1c, y1c, x2c ,y2c, SUB1=other_stars, MAX_DISTANCE=2.0
       if (other_stars[0] gt -1) then begin
        other_stars = image_model(x_near[other_stars],y_near[other_stars],f_near[other_stars],boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
;       if (debug gt 0) then begin
;        print, P
;        writefits, tmpdir + 'rawslice.fits', slice
;        writefits, tmpdir + 'other_stars.fits', image_shift(other_stars*P[1],P[2],P[3])
;        writefits, tmpdir + 'diff.fits', slice - image_shift(other_stars*P[1],P[2],P[3])
;        STOP
;       endif
        slice = slice - image_shift(P[0]+other_stars*P[1],P[2],P[3])
        endif
       endif
       stack[*,*,iref] = slice
      endfor
      w_frame = 0.0 ; initialize frame weight
      for iref = 0, nref-1 do begin
       subim = stack[*,*,iref]
       submask = psfmasks[*,*,iref]
       normim = subim - median(subim[bgind])
        ; the next IF statement avoids a problem
       ; if the entire subim = 0 (can occur when a reference
       ; star lies outside the crrent FOV
       if min(normim lt 0) then normim[where(normim lt 0)] = 0
       normim = circ_mask(normim, boxhw, boxhw, normrad,BORD=0)
       w_frame = w_frame + max(normim)
       if (unweighted eq 1 and total(normim) gt 0) then begin
        subim = subim/total(normim)
       endif
       if (subpix gt 0) then begin
        subim = image_shift(subim,xint[iref]-x_psf[iref],yint[iref]-y_psf[iref])
        submask = image_shift(submask,xint[iref]-x_psf[iref],yint[iref]-y_psf[iref])
        zeroind = where(submask lt 1, complement = ones)
        if (zeroind[0] gt -1) then submask[zeroind] = 0
        endif
       stack[*,*,iref] = subim
       psfmasks[*,*,iref] = submask
      endfor
      w_frame = 1.e-3 * w_frame/nref
;      print, 'Frame weight: ' + strn(w_frame)
      psf = stack_avg(stack, MASK=psfmasks, AVGYN=psfavg, CLIP=clip)

      ; noise cut and mask PSF
      ; ----------------------
      if (debug gt 0) then begin
       writefits, tmpdir + 'bgring_' +strtrim(string(it+1),2) + '.fits', psf*bgring
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'raw.fits', psf
      endif
;      print, 'Background: ' + string(median(psf[bgind]))
        noise = stddev(psf[bgind])
        psf = psf - median(psf[bgind])
        suppress = where(psf lt (nsigma[it+1] * noise), count)        
        if (count gt 0) then psf[suppress] = 0
        psfcen = centroid(circ_mask(psf,boxhw,boxhw,maskrad,BORDER=0))   
        psf = circ_mask(psf,psfcen[0],psfcen[1],maskrad,BORDER=0)
        psf = psf/total(psf)
        if (debug gt 0) then begin
         writefits, tmpdir + 'stack.fits', stack*psfmasks
         writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + '.fits', psf
       endif
     endfor ; end loop over iter
     temp = fltarr(nax1,nax2)
     temp[0:boxsize-1,0:boxsize-1] = psf
     psf = shift(temp,cen1-boxhw,cen2-boxhw)

  ; holography algorithm
  ; #####################

      ind = where(psf gt 0.0, count)

      if (count gt 0)  then begin

        if not(KEYWORD_SET(weightframes)) then w_frame = 1.0

        print, 'Frame weight: ' + strn(w_frame)

        ; image
        G = FFT(im)
        H = FFT(psf)*(nax1*nax2)
        HQ = CONJ(H)
        HABS = ABS(H)
        H2 = HABS^2
        GHQ = G*HQ

        hhh = REAL_PART(FFT(OTF*(GHQ/H2), /INVERSE))
        rcim = shift(hhh, center+1)

        if (debug gt 0) then begin
         writefits, tmpdir + 'rcim.fits', rcim
         STOP
        endif
       rcims[*,*,i] = rcim
 
       if keyword_set(psfout) then psfs[*,*,ngpsf] = psf
        ngpsf++

        if (withmask) then begin
         GM = FFT(thisfov)
         GMHQ = GM*HQ
         masknumer = masknumer + GMHQ * w_frame
         fimasknumer = fimasknumer + GMHQ * w_frame
        endif

        numer = numer + GHQ * w_frame
        denom = denom + H2 * w_frame
        nthis = nthis + 1
 
        finumer = finumer + GHQ * w_frame
        fidenom = fidenom + H2 * w_frame
        ngood++

        n123 = ngood mod nsub
        numersubreal[*,*,n123] = numersubreal[*,*,n123] + REAL_PART(GHQ) * w_frame
        denomsubreal[*,*,n123] = denomsubreal[*,*,n123] + REAL_PART(H2) * w_frame
        numersubim[*,*,n123] = numersubim[*,*,n123] + IMAGINARY(GHQ) * w_frame
        denomsubim[*,*,n123] = denomsubim[*,*,n123] + IMAGINARY(H2) * w_frame
        if (withmask) then begin
         masknumersubreal[*,*,n123] = masknumersubreal[*,*,n123] + REAL_PART(GMHQ) * w_frame
         masknumersubim[*,*,n123] = masknumersubim[*,*,n123] + IMAGINARY(GMHQ) * w_frame
        endif

        print, "Frame number: " + string(i+1) + " ok"

       endif else begin
         print, "Frame number: " + string(i+1) + " no speckle cloud detected"
       endelse   ; if (count gt 1.0)



    ; output of control data
    ; ----------------------

    if ((nthis mod out_iter) eq 0 and nthis gt 0) then begin
     hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
     rcim = shift(hhh, center+1)
     if (withmask) then begin
      mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
      rcmask = shift(mmm, center+1)
      support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
      rcim[support] = rcim[support]/rcmask[support]
      if (nosup[0] gt -1) then begin
       rcim[nosup] = 0
      endif
     endif
     writefits, tmpdir + 'rcim.fits', rcim
    endif
   endfor  ; end loop over frames in this cube
  

  ; output of PSFs
  if keyword_set(psfout) then writefits, '../psfs/psfs' + strtrim(string(ic+1),2) + '.fits', psfs[*,*,0:ngpsf-1]

  hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
  rcim = shift(hhh, center+1)
  if (withmask) then begin
   mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
   rcmask = shift(mmm, center+1)
   support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
   if (support[0] gt -1) then begin
    rcim[support] = rcim[support]/rcmask[support]
    if (nosup[0] gt -1) then begin
     rcim[nosup] = 0
     rcmask[nosup] = 0
    endif
   endif
  endif
  writefits, tmpdir + 'holo_' + strtrim(string(ic+1),2)+'.fits', rcim, header

  hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
  rcim = shift(hhh, center+1)
  if (withmask) then begin
   mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
   rcmask = shift(mmm, center+1)
   support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
   rcim[support] = rcim[support]/rcmask[support]
   if (nosup[0] gt -1) then begin
    rcim[nosup] = 0
    rcmask[nosup] = 0
   endif
   writefits, tmpdir + 'current_support.fits', rcmask
  endif
  writefits, tmpdir + 'current.fits', rcim

;  writefits, tmpdir + 'rcims' + strtrim(string(ic+1),2)+'.fits', rcims

endfor ; end loop over input cubes


; final image
; --------------------
hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
rcim = shift(hhh, center+1)
if (withmask) then begin
 mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
 rcmask = shift(mmm, center+1)
 support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
 rcim[support] = rcim[support]/rcmask[support]
 if (nosup[0] gt -1) then begin
  rcim[nosup] = 0.0
  rcmask[nosup] = 0
 endif
 writefits, outnam+'_support.fits', rcmask
endif
writefits, outnam+'.fits', rcim
print, 'Total number of frames used: ' + strtrim(string(ngood),2)

; sub-images
; --------------------

for is = 1, nsub do begin
 nhhh = complex(numersubreal[*,*,is-1],numersubim[*,*,is-1])
 dhhh = complex(denomsubreal[*,*,is-1],denomsubim[*,*,is-1])
 hhh = REAL_PART(FFT(OTF*(nhhh/dhhh), /INVERSE))
 rcim = shift(hhh, center+1)
 if (withmask) then begin
  maskhhh = complex(masknumersubreal[*,*,is-1],masknumersubim[*,*,is-1])
  mmm = REAL_PART(FFT(OTF*(maskhhh/dhhh), /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
  rcim[support] = rcim[support]/rcmask[support]
  if (nosup[0] gt -1) then begin
   rcim[nosup] = 0
   rcmask[nosup] = 0
  endif
 endif
 writefits, outnam+strtrim(string(is),2) + '.fits', rcim
endfor

print, 'Number of frames used in sub-images: ' + strtrim(string(float(ngood/nsub)),2)



END
 
