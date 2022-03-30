PRO CALEVAL, X, P, YMOD

  YMOD = P[0] + P[1]*X
  YMOD = image_shift(YMOD,P[2],P[3])

END


;---------------------

PRO EXTRACT_PSF, im, x_psf, y_psf, psfbox_hw, maskrad, normrad, NSIGMA=nsigma, REBFAC=rebfac, STARLIST=starlist, UNWEIGHTED = unweighted, DEBUG = debug, SUBPIX=subpix, ITER = iter

; VERSION DATE: 11 March 2014

; Define some defaut values
if  not keyword_set(nsigma) then nsigma = 3.0
if  not keyword_set(rebfac) then rebfac = 1.0
if  not keyword_set(unweighted) then unweighted = 0
if  not keyword_set(debug) then debug = 0
if  not keyword_set(subpix) then subpix = 1
if  not keyword_set(iter) then iter = 0

; number of iterations
 iter = n_elements(nsigma)

 ; star list
if keyword_set(starlist) then begin
  readcol, starlist, xx, yy, ff
  xx = rebfac * xx
  yy = rebfac * yy
  ff = rebfac^2 * ff
endif

; determine image size
sz = size(im)
nax1 = sz[1]
nax2 = sz[2]

; aplyh rebin factor
maskrad = rebfac*maskrad
nax1 = rebfac*nax1
nax2 = rebfac*nax2
normrad = rebfac*normrad
boxhw = rebfac*boxhw
im = FREBIN(im,nax1,nax2,/TOTAL)

; numbers for subarray extraction
boxsize = long(2*boxhw + 1)
cen1 = nax1/2
cen2 = nax2/2
center = [cen1,cen2]

; dummy needed  ringmask for background estimation
dummy = fltarr(boxsize,boxsize)
dummy[*,*] = 1
bgring = circ_mask(dummy,boxhw,boxhw,maskrad,BORDER=0, INNER=1)
bgind = where(bgring gt 0)

; reference sources
x_psf = rebfac * x_psf
y_psf = rebfac * y_psf
nref = n_elements(x_psf)
xint = round(x_psf)
yint = round(y_psf)

; PSF extraction
; ################

sub_arrays, im, xint, yint, boxsize, stack, masks
; save non-normalized stack for photometry
rawstack = stack

for iref = 0, nref-1 do begin
  subim = stack[*,*,iref]
  if (subpix gt 0) then begin
    subim = image_shift(subim,xint[iref]-x_psf[iref],yint[iref]-y_psf[iref])
  endif
  ; subtract local background and normalize sub-image 
  normim = subim - median(subim[bgind])
  normim = circ_mask(normim, boxhw, boxhw, normrad,BORD=0)
  if (unweighted eq 1 and total(normim) gt 0) then begin
    subim = subim/total(normim)
  endif
  stack[*,*,iref] = subim
 endfor
 psf = stack_avg(stack, MASK=masks, AVGYN=0)
 
 ; subtract remaining background, suppress noise,  and mask PSF
 ; -----------------------------------------------------------
   if (debug gt 0) then begin
     writefits, tmpdir + 'im.fits', im
     writefits, tmpdir + 'bgring_0.fits', psf*bgring
     writefits, tmpdir + 'psf_0raw.fits', psf
     writefits, tmpdir + 'stack_0.fits', stack*masks
    endif
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
  ; ----------------------------------

 for it = 0, iter-1 do begin

 ; clean surroundings of reference stars from secondary sources
 ; ------------------------------------------------------------
   stack = fltarr(boxsize,boxsize,nref)
   for iref = 0, nref-1 do begin
     sub_arrays, im, xint[iref], yint[iref], boxsize, slice, masks
     psfmasks[*,*,iref] = masks
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
      for iref = 0, nref-1 do begin
       subim = stack[*,*,iref]
       submask = psfmasks[*,*,iref]
       normim = subim - median(subim[bgind])
       normim = circ_mask(normim, boxhw, boxhw, normrad,BORD=0)
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
      psf = stack_avg(stack, MASK=psfmasks, AVGYN=psfavg, CLIP=clip)

      ; noise cut and mask PSF
      ; ----------------------
      if (debug gt 0) then begin
       writefits, tmpdir + 'bgring_' +strtrim(string(it+1),2) + '.fits', psf*bgring
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'raw.fits', psf
      endif
      noise = stddev(psf[bgind])
      psf = psf - median(psf[bgind]) - nsigma[it+1] * noise
      suppress = where(psf lt 0, count)        
      if (count gt 0) then psf[suppress] = 0
      psfcen = centroid(circ_mask(psf,boxhw,boxhw,maskrad,BORDER=0))   
      psf = circ_mask(psf,psfcen[0],psfcen[1],maskrad,BORDER=0)
      psf = psf/total(psf)
      if (debug gt 0) then begin
        writefits, tmpdir + 'stack.fits', stack*psfmasks
        writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + '.fits', psf
      endif
  endfor ; end loop over iter

END
 
