PRO CALEVAL, X, P, YMOD

  YMOD = P[0] + P[1]*X
  YMOD = image_shift(YMOD,P[2],P[3])

END


;---------------------

PRO WSA, indir, maskdir, masklist, inlist, outnam,  maskrad, normrad, nsigma, DEBUG = debug, rebfac, iter=iter, refsources, starlist, AIRY=airy, BORDER=bord, OUT_ITER=out_iter, PSFOUT = psfout, NORMPSF = normpsf, SUBPIX=subpix, PSFAVG=psfavg, CLIP=clip, tmpdir=tmpdir, BOXHW=boxhw


 ; star list
  readcol, starlist, xx, yy, ff

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

im = readfits(indir + list[0], NSLICE=0)
sz = size(im)
nax1 = sz[1]
nax2 = sz[2]

; multiply with rebin factor
maskrad = rebfac*maskrad
nax1 = rebfac*nax1
nax2 = rebfac*nax2
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
if (rebfac gt 1) then begin
  airy = CONGRID(airy,nax1,nax2,CUBIC=-0.5,/CENTER)
endif
airy = airy/total(airy)
APOD = FFT(airy)*(nax1*nax2)
OTF = ABS(APOD)

fimasknumer = fltarr(nax1,nax2)
finumer = fltarr(nax1,nax2)
fidenom = fltarr(nax1,nax2)
masknumer1 = fltarr(nax1,nax2)
numer1 = fltarr(nax1,nax2)
denom1 = fltarr(nax1,nax2)
masknumer2 = fltarr(nax1,nax2)
numer2 = fltarr(nax1,nax2)
denom2 = fltarr(nax1,nax2)
masknumer3 = fltarr(nax1,nax2)
numer3 = fltarr(nax1,nax2)
denom3 = fltarr(nax1,nax2)
ngood = 0L
n1 = 0L & n2 = 0L & n3 = 0L

; reference sources
x_psf = refsources[0,*]
y_psf = refsources[1,*]
f_psf = refsources[2,*]
nref = n_elements(x_psf)
xint = round(x_psf)
yint = round(y_psf)


for ic = 0, ncubes-1 do begin  ; start loop over all input cubes

  ; initialize variables for wsagraphy
  ;---------------------------------------

  masknumer = fltarr(nax1,nax2)
  numer = fltarr(nax1,nax2)
  numer[*,*] = 0.0
  denom = fltarr(nax1,nax2)
  denom[*,*] = 0.0
  nthis = 0L

  ; read speckle data cube
  ; --------------------------
  mcube = readfits(maskdir + mlist[ic])
  cube = readfits(indir + list[ic])
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
    thisfov = mcube[*,*,i]
  
    zeroind = where(thisfov lt 1, complement = ones)
    ; fov can have pixels with values << 1
    ; set them to zero to avoid problems with routines that
    ; use the masks
    if (zeroind[0] gt -1) then thisfov[zeroind] = 0
    thisfov[ones] = 1
    im = im * thisfov ; in case im has non-zero values outside of FOV
     
   ; rebin image and FOV mask
    if (rebfac gt 1.0) then begin
  ;  writefits, tmpdir + 'origim.fits', im
     im = CONGRID(im,nax1,nax2,/CENTER,CUBIC = -0.5)
     thisfov = CONGRID(thisfov,nax1,nax2,/CENTER,CUBIC = -0.5)
    endif
 
  ; PSF extraction
  ; ################

      psfim = im
      sub_arrays, psfim, xint, yint, boxsize, stack, masks
      ; save non-normalized stack for photometry
      rawstack = stack
      for iref = 0, nref-1 do begin
       subim = stack[*,*,iref]
       if (subpix gt 0) then begin
        subim = image_shift(subim,xint[iref]-x_psf[iref],yint[iref]-y_psf[iref])
       endif
       if (normpsf eq 1) then begin
        normim = circ_mask(subim, boxhw, boxhw, normrad,BORD=0)
        subim = subim/total(normim)
       endif
       stack[*,*,iref] = subim
      endfor
      psf = stack_avg(stack, MASK=masks, AVGYN=0)
 
       ; mask PSF
      ; -----------------------
      if (debug gt 0) then begin
       writefits, tmpdir + 'im.fits', im
       writefits, tmpdir + 'thisfov.fits', thisfov
       writefits, tmpdir + 'bgring_0.fits', psf*bgring
       writefits, tmpdir + 'psf_0raw.fits', psf
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
        writefits, tmpdir + 'stack_0.fits', stack
      endif

   
  ; improval of PSF with known sources
  ; ####################################

      for it = 0, iter-1 do begin

       ; clean surroundings of reference stars from secondary sources
       ; ------------------------------------------------------------
      stack = fltarr(boxsize,boxsize,nref)
      psfmasks = fltarr(boxsize,boxsize,nref)
      for iref = 0, nref-1 do begin
       sub_arrays, im, xint[iref], yint[iref], boxsize, slice, masks
       sub_arrays, thisfov, xint[iref], yint[iref], boxsize, maskslice, masks
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
        slice = slice - image_shift(other_stars*P[1],P[2],P[3])
        endif
       endif
       stack[*,*,iref] = slice
       psfmasks[*,*,iref] = maskslice
      endfor
      for iref = 0, nref-1 do begin
       subim = stack[*,*,iref]
       submask = psfmasks[*,*,iref]
       if (normpsf eq 1) then begin
        normim = circ_mask(subim, boxhw, boxhw, normrad,BORD=0)
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
;      print, 'Background: ' + string(median(psf[bgind]))
        noise = stddev(psf[bgind])
        psf = psf - (median(psf[bgind])+nsigma[it+1] * noise)  ; subtract background
        suppress = where(psf lt 0)        
        if (count gt 0) then psf[suppress] = 0
        psfcen = centroid(circ_mask(psf,boxhw,boxhw,maskrad,BORDER=0))   
;        psf = circ_mask(psf,psfcen[0],psfcen[1],maskrad,BORDER=0)
;        psf = psf/total(psf)
        if (debug gt 0) then begin
         writefits, tmpdir + 'stack.fits', stack
         writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + '.fits', psf
       endif
     endfor ; end loop over iter
     temp = fltarr(nax1,nax2)
     temp[0:boxsize-1,0:boxsize-1] = psf
     psf = shift(temp,cen1-boxhw,cen2-boxhw)

  ; wsagraphy algorithm
  ; #####################

      ind = where(psf gt 0.0, count)

      if (count gt 0)  then begin

        hmin = noise
        fwhm = fwhm(airy)
        sharplim = [0.01,1.0]
        roundlim = [-10.,10.]
        find, psf, xfound, yfound, ffound, sharp, roundness, hmin, fwhm, roundlim, sharplim
        print, xfound, yfound, ffound
        nfound = n_elements(ffound)
        psf[*,*] = 0.0
        xfound = round(xfound)
        yfound = round(yfound)
        for jf = 0, nfound-1 do begin
         psf[xfound[jf],yfound[jf]] = ffound[jf]
        endfor
        psf = psf/total(psf)
        if (debug gt 0) then begin
         writefits, tmpdir + 'speckles.fits', psf
        endif

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

        ; FOV mask
         GM = FFT(thisfov)
         GMHQ = GM*HQ

         masknumer = masknumer + GMHQ
         numer = numer + GHQ
         denom = denom + H2
         nthis = nthis + 1
         fimasknumer = fimasknumer + GMHQ
         finumer = finumer + GHQ
         fidenom = fidenom + H2
         ngood++

         n123 = ngood mod 3
         if (n123 eq 1) then begin
           masknumer1 = masknumer1 + GMHQ
           numer1 = numer1 + GHQ
           denom1 = denom1 + H2
           n1++
         endif else if (n123 eq 2) then begin
           masknumer2 = masknumer2 + GMHQ
           numer2 = numer2 + GHQ
           denom2 = denom2 + H2
           n2++
         endif else begin
          masknumer3 = masknumer3 + GMHQ
          numer3 = numer3 + GHQ
          denom3 = denom3 + H2
          n3++
         endelse

         print, "Frame number: " + string(i+1) + " ok"

       endif else begin
         print, "Frame number: " + string(i+1) + " no speckle cloud detected"
       endelse   ; if (count gt 1.0)



    ; output of control data
    ; ----------------------

    if ((nthis mod out_iter) eq 0 and nthis gt 0) then begin
     hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
     rcim = shift(hhh, center+1)
     mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
     rcmask = shift(mmm, center+1)
     support = where(abs(rcmask) gt 0.1, complement=nosup)
     rcim[support] = rcim[support]/rcmask[support]
     if (nosup[0] gt -1) then rcim[nosup] = 0.0
     writefits, tmpdir + 'rcim.fits', rcim
     if (nosup[0] gt -1) then rcim[nosup] = 0.0
;     writefits, tmpdir + 'support.fits', rcmask
    endif 
  endfor  ; end loop over frames in this cube
  

  ; output of PSFs
  if keyword_set(psfout) then writefits, '../psfs/psfs' + strtrim(string(ic+1),2) + '.fits', psfs[*,*,0:ngpsf-1]

  hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
  rcim = shift(hhh, center+1)
  mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(abs(rcmask) gt 0.1, complement=nosup)
  if (support[0] gt -1) then begin
   rcim[support] = rcim[support]/rcmask[support]
   if (nosup[0] gt -1) then rcim[nosup] = 0.0
   writefits, tmpdir + 'wsa_' + strtrim(string(ic+1),2)+'.fits', rcim
   if (nosup[0] gt -1) then rcim[nosup] = 0.0
;   writefits, tmpdir + 'support_' + strtrim(string(ic+1),2)+'.fits', rcmask
  endif
  hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
  rcim = shift(hhh, center+1)
  mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(abs(rcmask) gt 0.1, complement=nosup)
  writefits, tmpdir + 'current_support.fits', rcmask
  rcim[support] = rcim[support]/rcmask[support]
  if (nosup[0] gt -1) then rcim[nosup] = 0.0
  writefits, tmpdir + 'current.fits', rcim

;  writefits, tmpdir + 'rcims' + strtrim(string(ic+1),2)+'.fits', rcims

endfor ; end loop over input cubes


; final image
; --------------------
hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
rcim = shift(hhh, center+1)
mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
rcmask = shift(mmm, center+1)
support = where(abs(rcmask) gt 0.1, complement=nosup)
rcim[support] = rcim[support]/rcmask[support]
if (nosup[0] gt -1) then rcim[nosup] = 0.0
writefits, outnam+'.fits', rcim
if (nosup[0] gt -1) then rcim[nosup] = 0.0
writefits, outnam+'_support.fits', rcmask
print, 'Total number of frames used: ' + strtrim(string(ngood),2)

; sub-images
; --------------------
hhh = REAL_PART(FFT(OTF*(numer1/denom1), /INVERSE))
rcim = shift(hhh, center+1)
mmm = REAL_PART(FFT(OTF*(masknumer1/denom1), /INVERSE))
rcmask = shift(mmm, center+1)
support = where(abs(rcmask) gt 0.1, complement=nosup)
rcim[support] = rcim[support]/rcmask[support]
if (nosup[0] gt -1) then rcim[nosup] = 0.0
writefits, outnam+'1.fits', rcim
print, 'Number of frames used in sub-image 1: ' + strtrim(string(n1),2)
hhh = REAL_PART(FFT(OTF*(numer2/denom2), /INVERSE))
rcim = shift(hhh, center+1)
mmm = REAL_PART(FFT(OTF*(masknumer2/denom2), /INVERSE))
rcmask = shift(mmm, center+1)
support = where(abs(rcmask) gt 0.1, complement=nosup)
rcim[support] = rcim[support]/rcmask[support]
if (nosup[0] gt -1) then rcim[nosup] = 0.0
writefits, outnam+'2.fits', rcim
print, 'Number of frames used in sub-image 2: ' + strtrim(string(n2),2)
hhh = REAL_PART(FFT(OTF*(numer3/denom3), /INVERSE))
rcim = shift(hhh, center+1)
mmm = REAL_PART(FFT(OTF*(masknumer3/denom3), /INVERSE))
rcmask = shift(mmm, center+1)
support = where(abs(rcmask) gt 0.1, complement=nosup)
rcim[support] = rcim[support]/rcmask[support]
if (nosup[0] gt -1) then rcim[nosup] = 0.0
writefits, outnam+'3.fits', rcim
print, 'Number of frames used in sub-image 3: ' + strtrim(string(n3),2)


END
