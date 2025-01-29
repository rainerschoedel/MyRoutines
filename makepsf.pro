PRO CALEVAL, X, P, YMOD

  YMOD = P[0]*X
  YMOD = image_shift(YMOD,P[1],P[2])

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
;            The default value of MINDIST is 2 pixels.
; LAST MODIFICATION
;
; Introduced MINDIST keyword to avoid that stars very close to the
; reference star (or the reference star itself) can be subtracted as secondaries.
; 
; Activated the CALEVAL procedure again that scales the photometry of
; secondaries to the current best-estimate of the PSF.
; Rainer Schoedel 4 March 2015
;
; Changed nme from MYPSF_REPEAT_EXTRACT to MAKEPSF
; Changed normalisation and weighting of sources.
; Can also be used for first estimation of PSF
; Removed determination of constant background from CALEVAL 
; because it can create problems with nearby unresolved 
; sources or nearby sources that are not included in the list of
; secondaries
; Included potential use of mask
; PSF must be fully covered (important in case of unregularly masked
; image
; Note that iterations appear to make the result slightly worse
;
; 24 July 2024, Rainer Schoedel
; USE_CENTROID can be set to sub-pixel shift the reference stars with the 
; StarFinder CENTROID function. This is necessary if the positions
; of the reference stars have not yet been determined with high accuracy.
; 
;

PRO MAKEPSF, x_ref, y_ref, x_stars, y_stars, f_stars, image, nrad, FOVMASK = thisfov, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_sigma, MASKRAD = mrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=local_sky, USE_CENTROID=use_centroid


;  print, 'HERE'
 if not(keyword_set(debug)) then debug = 0
 if not(keyword_set(mindist)) then mindist = 2.
 if keyword_set(mask) then withmask = 1 else withmask = 0
 if (not KEYWORD_SET(satlevel)) then satlevel = 1.e9
 if (not KEYWORD_SET(local_sky)) then local_sky = 0
 if (not KEYWORD_SET(thisfov)) then begin
    thisfov = image
    thisfov[*,*] = 1
 endif
 if not(keyword_set(use_centroid)) then use_centroid = 0

 sz = size(psf)
 boxsize = sz[1]
 boxhw = boxsize/2 ; careful: MUST NOT BE FLOAT!
 
 sz = size(image)
 nax1 = sz[1]
 nax2= sz[2]

 nref = n_elements(x_ref)
 xint = round(x_ref)
 yint = round(y_ref)

 psfim = image

if KEYWORD_SET(background) then psfim = psfim - background

bgring = replicate(1,boxsize,boxsize)
bgring = CIRC_MASK(bgring,boxhw,boxhw,mrad,/INNER)
bgind = where(bgring gt 0, n_bgind)

dummy = replicate(1,boxsize,boxsize)
dummy = CIRC_MASK(bgring,boxhw,boxhw,mrad)
refind = where(dummy gt 0, npix_ref)

; 1) Select stars in FoV
;use psfmasks if a reference star may not be contained in all
; cubes because of dithering
; ---------------------------

  sub_arrays, psfim, xint, yint, boxsize, stack, masks
  if withmask then begin
    sub_arrays, thisfov, xint, yint, boxsize, psfmasks, masks
  endif
  if debug then begin
    writefits, tmpdir + 'rawstack.fits', stack
;    writefits, 'masks.fits', masks
  endif
  valid_ref = replicate(1,nref)
  for iref = 0, nref-1 do begin
   if withmask then masks[*,*,iref] = masks[*,*,iref]*psfmasks[*,*,iref]
   dummy = circ_mask(masks[*,*,iref],boxhw,boxhw,mrad)
   if total(dummy) lt (0.99*npix_ref) then begin
     valid_ref[iref] = 0
   endif
  endfor
  acceptref = where(valid_ref gt 0, nref_accept)
  xint_accept = xint[acceptref]
  yint_accept = yint[acceptref]
  x_psf_accept = x_ref[acceptref]
  y_psf_accept = y_ref[acceptref]
  stack = stack[*,*,acceptref]
  masks = masks[*,*,acceptref]

  print, 'Using: ' + strn(nref_accept) + ' reference stars. '
  psf_weights = fltarr(nref_accept)
  ; save current stack in rawstack
  rawstack = stack

; (1) Use current PSF and information on image
; to improve PSF estimate
; -----------------------------------------
if (iter gt 0) then begin

  psf = psf/total(psf)

 ; improval of PSF with known sources
  ; ####################################

  for it = 0, iter-1 do begin

    ; clean surroundings of reference stars from secondary sources
    ; ------------------------------------------------------------

    for iref = 0, nref_accept-1 do begin
      slice = rawstack[*,*,iref]
      ; estimate local sky
      if local_sky then begin
        mmm, slice[bgind], skymod, skysig, skyskew, /SILENT 
        slice = slice - skymod     
      endif
 ;       print, iref
 ;       print, skymod
      distances = sqrt((x_stars-x_psf_accept[iref])^2 + (y_stars-y_psf_accept[iref])^2)
      nearby = where(distances lt sqrt(2.)*boxsize/2.)
      if (nearby[0] gt -1) then begin
        x_near = x_stars[nearby] - (xint_accept[iref] - boxhw)
        y_near = y_stars[nearby] - (yint_accept[iref] - boxhw)
        f_near = f_stars[nearby]
        ; sort by decreasing flux
        ; needed for latter residual masking of brightest stars
        ord = reverse(sort(f_near))
        x_near = x_near[ord]
        y_near = y_near[ord]
        f_near = f_near[ord]
        model = image_model(x_near,y_near,f_near,boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
        P = [1.,0.,0.]
        W = slice
        W = sqrt(abs(W))
;        W[*,*] = 1.0
        W = W * masks[*,*,iref]
        res = mpcurvefit(model,slice,W,P,sigma,FUNCTION_NAME='CALEVAL',/NODERIVATIVE,/QUIET)
        compare_lists, x_near, y_near, x_psf_accept[iref]-(xint_accept[iref]-boxhw), y_psf_accept[iref]-(yint_accept[iref]-boxhw), x1c, y1c, x2c ,y2c, SUB1=other_stars, MAX_DISTANCE=mindist
       if (other_stars[0] gt -1) then begin
         secondaries = image_model(x_near[other_stars],y_near[other_stars],f_near[other_stars],boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
          
;           if (debug gt 0) then begin
;             if local_sky then print, skymod
;             print, P
;             print, iref
;             writefits, tmpdir + 'model.fits', model
;             writefits, tmpdir + 'W.fits', W
;             writefits, tmpdir + 'rawslice.fits', slice
;             writefits, tmpdir + 'secondaries.fits', image_shift(secondaries*P[0],P[1],P[2])
;             writefits, tmpdir + 'diff.fits', slice - image_shift(secondaries*P[0],P[1],P[2])
;             STOP
;            endif
          slice = slice - image_shift(secondaries*P[0],P[1],P[2])
          endif       ;  (other_stars[0] gt -1)
        endif        ;  if (nearby[0] gt -1)
        stack[*,*,iref] = slice
       endfor

       ; Align with sub-pixel shift
       for iref = 0, nref_accept-1 do begin
        subim = stack[*,*,iref]
        submask = masks[*,*,iref]
        subim = image_shift(subim,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
        submask = image_shift(submask,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
        zeroind = where(submask lt 0.99, complement = ones)
        if (use_centroid) then begin
          subim = CENTROIDER(subim,XSHIFT=xs,YSHIFT=ys)
          submask = image_shift(submask,xs,ys)
        endif
        if (zeroind[0] gt -1) then begin
          submask[zeroind] = 0
          submask[ones] = 1
        endif  
        saturated = where(subim gt satlevel)
        submask[saturated] = 0
        normim = subim
        if min(normim lt 0) then normim[where(normim lt 0)] = 0
        normim = circ_mask(normim, boxhw, boxhw, nrad)
        ; normalization will underweight stars with saturated cores
        ; should be no problem
         fnorm = total(normim)
         subim = subim/fnorm * submask
         psf_weights[iref] = sqrt(fnorm)
         stack[*,*,iref] = subim
         masks[*,*,iref] = submask
       endfor
       psf_weights = round(psf_weights/min(psf_weights))
       if (unweighted eq 1) then psf_weights[*] = 1
       print, 'Weights: '
       print, psf_weights
       psf = stack_median(stack, WEIGHTS=psf_weights, MASK=masks)
       psf_sigma = stack_error(stack, WEIGHTS=psf_weights, MASK=masks)
       
      ; suppress noisy parts and mask PSF, obtain mean S/N
      ; ---------------------------------------------------
      if (debug gt 0) then begin
       writefits, tmpdir + 'bgring_' +strtrim(string(it+1),2) + '.fits', psf*bgring
       writefits, tmpdir + 'stack.fits', stack
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'raw.fits', psf
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'sigma.fits', psf_sigma
      endif


 ;      if debug then writefits, tmpdir + 'psf_nomask.fits', psf
 ;      print, 'iter:  ' + strn(iter)
 ;      MASK_PSF, psf, mrad, PSF_MASKED=psf_masked ; to avoid having negative wings of the PSF
 ;      psf = psf_masked
 ;      if debug then writefits, tmpdir + 'psf_masked.fits', psf
 ;      psfnorm = total(psf)
 ;      psf = psf/psfnorm

    endfor                      ; end loop over iter

; if PSF does not yet exist, extract one
; ------------------------------------------
endif else begin

  psf_weights = intarr(nref_accept) ; for optional weighting
  for iref = 0, nref_accept-1 do begin
    subim = stack[*,*,iref]
    submask = masks[*,*,iref]
    ; estimate local sky
    if local_sky then begin
      mmm, subim[bgind], skymod, skysig, skyskew, /SILENT       
      subim = subim - skymod
    endif
    ; sub-pixel shift
    subim = image_shift(subim,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
    submask = image_shift(submask,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
    zeroind = where(submask lt 0.99, complement = ones)
    if (use_centroid) then begin
      subim = CENTROIDER(subim,XSHIFT=xs,YSHIFT=ys)
      submask = image_shift(submask,xs,ys)
    endif
    if (zeroind[0] gt -1) then begin
      submask[zeroind] = 0
      submask[ones] = 1
    endif  
    subim = subim * submask
    masks[*,*,iref] = submask

    ; Normalize the reference stars
    ; normalization will underweight stars with saturated cores
    ; should be no problem
    normim = subim
    normim = circ_mask(normim, boxhw, boxhw, nrad)
    fnorm = total(normim)
    saturated = where(subim gt satlevel)
    submask[saturated] = 0
    normim = normim * submask
    subim = subim/fnorm
    psf_weights[iref] = sqrt(fnorm)
    stack[*,*,iref] = subim*submask
    masks[*,*,iref] = submask
   endfor
   psf_weights = round(psf_weights/min(psf_weights))
   if (unweighted eq 1) then psf_weights[*] = 1
   print, 'Weights: '+ strn(psf_weights)
   psf = stack_median(stack, WEIGHTS=psf_weights, MASK=masks)
   psf_sigma = stack_error(stack, WEIGHTS=psf_weights, MASK=masks)

   if (debug gt 0) then begin
     writefits, tmpdir + 'im.fits', psfim
     if withmask then writefits, tmpdir + 'thisfov.fits', thisfov
     writefits, tmpdir + 'bgring.fits', psf*bgring
     writefits, tmpdir + 'psfraw.fits', psf
     writefits, tmpdir + 'rawstack.fits', rawstack
     writefits, tmpdir + 'stackmasks.fits', masks
   endif

   ; circular mask and normalise PSF
   ; -------------------------------
;   MASK_PSF, psf, mrad, PSF_MASKED=psf_masked ; to avoid having negative wings of the PSF
;   psf = psf_masked
;   psf = psf/total(psf)
 
endelse


END
