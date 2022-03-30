PRO CALEVAL, X, P, YMOD

  YMOD = P[0]*X
  YMOD = image_shift(YMOD,P[1],P[2])

END


;---------------------

PRO HOLO_MOSAIC, indir, inlist, outnam,  maskrad, normrad, nsigma, rebfac, refsources, starlist, maskdir=maskdir, masklist=masklist, DEBUG = debug, iter=iter, AIRY=airy, OUT_ITER=out_iter, PSFOUT = psfout, UNWEIGHTED = unweighted, SUBPIX=subpix, tmpdir=tmpdir, BOXHW=boxhw, NSUB=nsub, MINSUPP=minsupp, MAXSUPP = maxsupp,   REBITER=rebiter, RAWOUT = rawout, N_MASK_SECONDARY=n_mask_secondary, PSF_FRAC = psf_frac, CORRECT_SKY = correct_sky, SMOOTHMASK = smoothmask, N_REF_MAX = n_ref_max, PR = pr, CIRC_BORDER=circ_border, SATLEVEL=satlevel, PSF_BORDER=psf_border, ESTIM_BG = estim_bg, SAT_MASK = sat_mask

; ---------------------------------------------------
; As HOLO_FOURIER, but subtracts the noise threshold
; ---------------------------------------------------

withmask = KEYWORD_SET(masklist)
if (not KEYWORD_SET(debug)) then debug = 0
if (not KEYWORD_SET(maxsupp)) then maxsupp = 1.1
if (not KEYWORD_SET(minsupp)) then minsupp = 0.9
if (not KEYWORD_SET(n_mask_secondary)) then n_mask_secondary = 0
if (not KEYWORD_SET(rawout)) then rawout= 0
if (not KEYWORD_SET(rebiter)) then rebiter = 0
if (not KEYWORD_SET(psf_frac)) then psf_frac = 0.9
if (not KEYWORD_SET(correct_sky)) then correct_sky = 0
if (not KEYWORD_SET(smoothmask)) then smoothmask = 0
if (not KEYWORD_SET(pr)) then pr = 1
if (not KEYWORD_SET(circ_border)) then circ_border = 0
if (not KEYWORD_SET(satlevel)) then satlevel = 1.0e9
if (not KEYWORD_SET(psf_border)) then psf_border = 0.0
if (not KEYWORD_SET(estim_bg)) then estim_bg = 0


; Identical to holo7sub.pro except that it 
; subtracts the noise threshold from the final PSF
; -------------------------------------------------

 ; star list
  readcol, starlist, xx, yy, ff
 ; rebin positions if necessary
 if (rebiter lt 1) then begin
  xx = rebfac * xx
  yy = rebfac * yy
  ff = rebfac^2 * ff
 endif
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
satlev = satlevel/rebfac^2 ; variable must be changed
                                ; otherwise satlevel is passed back
                                ; with changed value, which will 
                                ; create problems in repeated calls of
                                ; this routine.
maskrad = rebfac * maskrad
nax1 = rebfac*nax1
nax2 = rebfac*nax2
normrad = rebfac*normrad
boxhw = rebfac*boxhw
; pixel indixes for subarray extraction
boxsize = long(2*boxhw + 1)

; number of pixels in image
;npix = float(nax1 * nax2)

cen1 = nax1/2
cen2 = nax2/2
center = [cen1,cen2]

dmax = rebfac * 2.0  ; contaminating  stars close to a reference stars will be
                     ; recognized if they are located at more than dmax 
                     ; from the reference star
                     ; smaller dmax values can be dangerous because
                     ; then a reference star may be considered a contaminator
                     ; which would screw up the PSF extraction

; dummy needed  ringmask for background estimation
; npix_ref if the minimum required support region for
; reference stars (in case of variable masking/dithering)
dummy = fltarr(boxsize,boxsize)
dummy[*,*] = 1
bgring = circ_mask(dummy,boxhw,boxhw,maskrad, INNER=1)
bgind = where(bgring gt 0, n_bgind)

; SUPPORT FOR PSF, NOT IN PYTHON CODE
dummy[*,*] = 1
dummy = circ_mask(dummy,boxhw,boxhw,maskrad)
psf_support = where(dummy gt 0, npix_ref)

; ######################################

; apodization with telescope OTF
; be careful to provide an adequately sampled PSF 
; of the correct size 
; (re-binning of PSF can be done but is non-optimal)
airy = airy/total(airy)
APOD = FFT(airy)*(nax1*nax2)
OTF = ABS(APOD)

nexp = fltarr(nax1,nax2) ; records number of exposures per pixel
finumer = fltarr(nax1,nax2)
fidenom = fltarr(nax1,nax2)
finumersubreal = fltarr(nax1,nax2,nsub)
fidenomsubreal = fltarr(nax1,nax2,nsub)
finumersubim = fltarr(nax1,nax2,nsub)
fidenomsubim = fltarr(nax1,nax2,nsub)
ngood = 0L

fimasknumer = fltarr(nax1,nax2)
fimasknumersubreal = fltarr(nax1,nax2,nsub)
fimasknumersubim = fltarr(nax1,nax2,nsub)

; reference sources
x_psf = refsources[0,*]
y_psf = refsources[1,*]
f_psf = refsources[2,*]
; rebin positions if neessary
if (rebiter lt 1) then begin
 x_psf = rebfac * x_psf
 y_psf = rebfac * y_psf
 f_psf = rebfac^2 * f_psf
endif
nref = n_elements(x_psf)
xint = round(x_psf)
yint = round(y_psf)

; open files for output lists
openw, outfiles1, tmpdir + 'holo_ims.txt', /get_lun
openw, outfiles2, tmpdir + 'im_masks.txt', /get_lun
openw, outfiles3, tmpdir + 'weights.txt', /get_lun

; delta function
delta_width = 7*pr
delta = fltarr(delta_width,delta_width)
delta_cen = delta_width/2
delta[delta_cen,delta_cen] = 1.0



; Now start lop over all input cubes
; --------------------------------
for ic = 0, ncubes-1 do begin 

  ; initialize variables for holography
  ;---------------------------------------

  expmap = replicate(0.0,nax1,nax2) ; records number of exposures per pixel
  mask = replicate(0.0,nax1,nax2)   ; records edge effects of holography
  masknumer = replicate(0.0,nax1,nax2)
  numer = replicate(0.0,nax1,nax2)
  denom = replicate(0.0,nax1,nax2)
  numersubreal = fltarr(nax1,nax2,nsub)
  denomsubreal = fltarr(nax1,nax2,nsub)
  numersubim = fltarr(nax1,nax2,nsub)
  denomsubim = fltarr(nax1,nax2,nsub)
  masknumersubreal = fltarr(nax1,nax2,nsub)
  masknumersubim = fltarr(nax1,nax2,nsub)
  nthis = 0L

  ; read speckle data cube
  ; --------------------------
  if (withmask) then begin
    mcube = readfits(maskdir + mlist[ic])
    thisfov = mcube[*,*,0]  
     if (rebfac gt 1) then  thisfov = CREBIN(thisfov,nax1,nax2) ; Do not normalize rebinned FOV!
     ; next two lines are necessary
     ; if the mask has previously been shifted by sub-pixels
     ; or if rebinning has created pixels with values < 1
     edgpix = where(thisfov ne 1)
     if edgpix[0] gt -1 then thisfov[edgpix] = 0
     if SMOOTHMASK gt 0 then begin
       smoothfov = gauss_smooth(thisfov,smoothmask,/EDGE_TRUNCATE)
       smoothfov = smoothfov * thisfov
       thisfov = smoothfov
       ; repeat smoothing to obtain really
       ; smooth edges (otherwise there is a jump near the edge)
       smoothfov = gauss_smooth(thisfov,smoothmask,/EDGE_TRUNCATE)
       smoothfov = smoothfov * thisfov
       thisfov = smoothfov
       smoothfov = gauss_smooth(thisfov,smoothmask,/EDGE_TRUNCATE)
       smoothfov = smoothfov * thisfov
       thisfov = smoothfov
     endif
   endif else thisfov = replicate(1,nax1,nax2) ; define thisfov for exposure map

   cube = readfits(indir + list[ic], header)
   sz = size(cube)
   if (sz[0] gt 2) then  nim = sz[3] else nim = 1
   print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
   print, 'Read in cube: ' + list[ic]
   print, 'Number of frames in this cube: ' + strtrim(string(nim),2)
   print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'

  ; Select reference stars in FoV
  ; mask some areas
  ; create delta map of reference stars for
  ; initial PSF estimation
  ; --------------------------------------

  ; Valid reference sources
  ; are those sources that are contained
  ; in field of view (at least to psf_frac part)
  valid_ref = replicate(1,nref)
  sub_arrays, thisfov, xint, yint, boxsize, psffovs, masks
  for iref = 0, nref-1 do begin
    masks[*,*,iref] = masks[*,*,iref]*psffovs[*,*,iref]
    dummy = circ_mask(masks[*,*,iref],boxhw,boxhw,maskrad)
    if total(dummy) lt round(psf_frac*npix_ref) then begin
      valid_ref[iref] = 0
    endif
  endfor
  acceptref = where(valid_ref gt 0, nref_accept,complement=reject)

  if nref_accept gt 0 then begin
    xint_accept = xint[acceptref]
    yint_accept = yint[acceptref]
    x_psf_accept = x_psf[acceptref]
    y_psf_accept = y_psf[acceptref]
    f_psf_accept = f_psf[acceptref]
    xint_reject = xint[reject]
    yint_reject = yint[reject]

    if KEYWORD_SET(n_ref_max) and (n_ref_max lt nref_accept) then begin
     xint_accept = xint_accept[0:n_ref_max-1]
     yint_accept = yint_accept[0:n_ref_max-1]
     x_psf_accept = x_psf_accept[0:n_ref_max-1]
     y_psf_accept = y_psf_accept[0:n_ref_max-1]
     f_psf_accept = f_psf_accept[0:n_ref_max-1]
     nref_accept = n_ref_max
    endif
    print, 'Using: ' + strn(nref_accept) + ' reference stars. '

 
    if keyword_set(psfout) then psfs = fltarr(nax1,nax2,nim)
    ngpsf = 0L

  ; start loop over frames in this cube
  ; ----------------------------------
  for i = 0L, nim-1 do begin

    ; load ith image 
    im = cube[*,*,i]
    ; Rebin image 
    if (rebfac gt 1) then im = CREBIN(im,nax1,nax2,/TOTAL)
    im = im * thisfov
    expmap = expmap + thisfov

 
  ; PSF extraction
  ; ################
     psfim = im

     ; First PSF extraction
     ; Will be unweighted.
     ; ------------------

     sub_arrays, psfim, xint_accept, yint_accept, boxsize, stack, masks
     sub_arrays, thisfov, xint_accept, yint_accept, boxsize, psfmasks, masks
     masks = masks * psfmasks

     for iref = 0, nref_accept-1 do begin
      subim = stack[*,*,iref]
      submask = masks[*,*,iref]
      sky_pixels = subim[bgind]
      sky_mask = submask[bgind]
      sky_good = where(sky_mask gt 0)
      bg_level = median(sky_pixels[sky_good])
      normim = subim - bg_level
      if (subpix gt 0) then begin
       subim = image_shift(subim,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
      endif
      if min(normim lt 0) then normim[where(normim lt 0)] = 0
      normim = circ_mask(normim, boxhw, boxhw, normrad)
      subim = subim/total(normim)
      stack[*,*,iref] = subim
     endfor
     psf = stack_median(stack, MASK=masks)
     if (debug gt 0) then writefits, tmpdir + 'psf_0raw.fits', psf

     psf_sigma = psf
     RESISTANT_Mean,psf[bgind],3.0,mean_psf,sigma_psf,Num_Rej
     noise = sigma_psf * sqrt(n_bgind - Num_Rej)
     psf = psf - mean_psf- nsigma[0] * noise
     suppress = where(psf lt 0, count)        
     if (count gt 0) then psf[suppress] = 0
;     psfcen = centroid(circ_mask(psf,boxhw,boxhw,,maskrad))
     psf = circ_mask(psf,boxhw,boxhw,maskrad)
     psf = psf/total(psf)
     ; compute PSF FWHM for later use
     psf_fwhm = round(fwhm(psf))
     psf_sigma[*,*] = noise

     if (debug gt 0) then begin
       writefits, tmpdir + 'im.fits', im
       writefits, tmpdir + 'thisfov.fits', thisfov
       writefits, tmpdir + 'stack0.fits', stack
       writefits, tmpdir + 'bgring_0.fits', psf*bgring
       writefits, tmpdir + 'psf_0.fits', psf
    endif
  
  ; improval of PSF with known sources
  ; ####################################

 
     psf_weights = intarr(nref_accept) ; for optional weighting

      ; iterative improval of PSF
      for it = 0, iter-1 do begin

       ; clean surroundings of reference stars from secondary sources
       ; ------------------------------------------------------------
      stack = fltarr(boxsize,boxsize,nref_accept)
      psfmasks = fltarr(boxsize,boxsize,nref_accept)

     ; Estimate and subtract diffuse background
     ; ----------------------------------------

      model = image_model(xx,yy,ff,nax1,nax2,psf,REFERENCE_PIX=[boxhw,boxhw],FTOL=1.D-2)
      model = model * thisfov
      P = [1.,0.,0.]
      W = model
      gt0 = where(model gt 0,complement=lt0)
      W[gt0] = sqrt(W[gt0])
      W[lt0] = 0
      res = mpcurvefit(model,im,W,P,sigma,FUNCTION_NAME='CALEVAL',/NODERIVATIVE,/QUIET)
;      print, 'Model fitting parameters for overall image: ' + strn(P)
      im_model = image_shift(model*P[0],P[1],P[2])
      diffuse = im - im_model
      if estim_bg then begin
         background = estimate_background(diffuse, maskrad, CUBIC=-0.5)
         ; estimate overall sky (if this frame has an offset)
          mmm, background, skymod, skysig, skyskew, /SILENT       
      endif else begin
          mmm, im, skymod, skysig, skyskew, /SILENT        
          background = im
          background[*,*] = skymod
      endelse

      if (debug gt 0) then begin
        print, P
        writefits, tmpdir + 'model.fits', im_model
        writefits, tmpdir + 'W.fits', W
      endif

      ; use psfmasks if a reference star may not be contained in all
      ; cubes because of dithering
      ; ---------------------------

      for iref = 0, nref_accept - 1 do begin
       sub_arrays, background, xint_accept[iref], yint_accept[iref], boxsize, thisbg, masks
       sub_arrays, im, xint_accept[iref], yint_accept[iref], boxsize, slice, masks
       slice = slice - thisbg
       sub_arrays, thisfov, xint_accept[iref], yint_accept[iref], boxsize, maskslice, masks
       psfmasks[*,*,iref] = maskslice * masks
       distances = sqrt((xx-x_psf_accept[iref])^2 + (yy-y_psf_accept[iref])^2)
       nearby = where(distances lt sqrt(2.)*boxsize/2.)
       if (nearby[0] gt -1) then begin
         x_near = xx[nearby] - (xint_accept[iref] - boxhw)
         y_near = yy[nearby] - (yint_accept[iref] - boxhw)
         f_near = ff[nearby]
         ; sort by decreasing flux
         ; needed for latter residual masking of brightest stars
         ord = reverse(sort(f_near))
         x_near = x_near[ord]
         y_near = y_near[ord]
         f_near = f_near[ord]
         compare_lists, x_near, y_near, x_psf_accept[iref]-(xint_accept[iref]-boxhw), y_psf_accept[iref]-(yint_accept[iref]-boxhw), x1c, y1c, x2c ,y2c, SUB1=other_stars, MAX_DISTANCE=dmax 
         if (other_stars[0] gt -1) then begin
           secondaries = image_model(x_near[other_stars],y_near[other_stars],f_near[other_stars],boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
          
           slice = slice - (background + image_shift(secondaries*P[0],P[1],P[2]))
;         print, 'P[0]: ' + strn(P[0])
         ; subtraction of secondary sources will not be perfect
         ; to further suppress there influence
         ; mask small circular regions centred on the secondaries
         ; mask only the remnants of the brightest secondaries
         ; otherwise this becomes a problem in crowded fields
         if (n_mask_secondary gt 0) then begin
           n_mask = min([n_elements(other_stars),n_mask_secondary])
            ; The PSF will generally not be centered on the center of 
            ; array. Therefore we have to do some acrobatics here.
           psfmax = max(psf,max_sub)
           xymax = array_indices(psf,max_sub)
           psf_mask = circ_mask(psf,xymax[0],xymax[1],round(psf_FWHM/2.))
           x_mask = round((x_near[other_stars[0:n_mask-1]]))
           y_mask = round((y_near[other_stars[0:n_mask-1]]))
           f_mask = round((f_near[other_stars[0:n_mask-1]]))  ; just needed for formal reasons
           mask_secondaries = image_model(x_mask,y_mask,f_mask,boxsize,boxsize,psf_mask,REFERENCE_PIX=[boxhw,boxhw])
           mask_ind = where(mask_secondaries gt 0, complement=nomask)
           maskslice = psfmasks[*,*,iref]
           if (mask_ind[0] gt -1) then maskslice[mask_ind] = 0
             psfmasks[*,*,iref] = maskslice
           endif      ;  if (n_mask_secondary gt 0)
          endif           ;  (other_stars[0] gt -1)
;          if debug then begin
;             print, 'Number of nearby stars: ' + strn(n_elements(nearby))
;             writefits, tmpdir + 'slice.fits', slice
;             writefits, tmpdir + 'secondaries.fits', secondaries
;             STOP
;          endif
        endif        ;  if (nearby[0] gt -1)
        stack[*,*,iref] = slice
       endfor
      
       if (debug gt 0) then writefits,  tmpdir + 'rawstack_1.fits', stack
       for iref = 0, nref_accept-1 do begin
         subim = stack[*,*,iref]
         submask = psfmasks[*,*,iref]
         if (subpix gt 0) then begin
          subim = image_shift(subim,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
          submask = image_shift(submask,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
          zeroind = where(submask lt 0.99, complement = ones)
          if (zeroind[0] gt -1) then begin
           submask[zeroind] = 0
           submask[ones] = 1
          endif  
         endif
         saturated = where(subim gt satlev)
         submask[saturated] = 0
         normim = subim
         if min(normim lt 0) then normim[where(normim lt 0)] = 0
         normim = circ_mask(normim, boxhw, boxhw, normrad)
        ; normalization will underweight stars with saturated cores
        ; should be no problem
         fnorm = total(normim)
         subim = subim/fnorm * submask
         psf_weights[iref] = sqrt(fnorm)
         stack[*,*,iref] = subim
         psfmasks[*,*,iref] = submask
       endfor
       psf_weights = round(psf_weights/min(psf_weights))
       if (unweighted eq 1) then psf_weights[*] = 1
       psf = stack_median(stack, WEIGHTS=psf_weights, MASK=psfmasks)
       psf_sigma = stack_error(stack, MASK=psfmasks)

      ; suppress noisy parts and mask PSF, obtain mean S/N
      ; ---------------------------------------------------
      if (debug gt 0) then begin
       writefits, tmpdir + 'bgring_' +strtrim(string(it+1),2) + '.fits', psf*bgring
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'raw.fits', psf
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'sigma.fits', psf_sigma
      endif
;      print, 'Background: ' + string(median(psf[bgind]))
;      use astrolib robus statistics script
;      to estimate mean and stddev
;      this helps to avoid rare problems with individual bright pixels
       RESISTANT_Mean,psf[bgind],3.0,mean_bg,sigma_bg,Num_Rej
       noise = sigma_bg * sqrt(n_bgind - Num_Rej)
       psf = psf - mean_bg - nsigma[it+1] * noise
       suppress = where(psf  lt 0, complement=accept,count)
       if (count gt 0) then psf[suppress] = 0
       psf = circ_mask(psf,boxhw,boxhw,maskrad)
       psfnorm = total(psf)
       psf = psf/psfnorm
       if (debug gt 0) then begin
        writefits, tmpdir + 'stack.fits', stack*psfmasks
        writefits, tmpdir + 'stackmasks.fits', psfmasks
        writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + '.fits', psf
     endif
     endfor ; end loop over iter
     ; adapt PSF to size of input frame
     ; shift it to the center (only integer shift!)
     temp = fltarr(nax1,nax2)
     temp[0:boxsize-1,0:boxsize-1] = psf
     psf = shift(temp,cen1-boxhw,cen2-boxhw)
 
; holography algorithm
  ; #####################

      ind = where(psf gt 0.0, count)

      if (count gt 0)  then begin
         
        if (correct_sky gt 0) then begin
         print, 'sky from PSF: ' + strn(sky)
         im = thisfov*(im - skymod)  ; subtract sky
        endif
        im = im * thisfov  ; for correct weighting of mask edges

        G = FFT(im)
        H = FFT(psf)*(nax1*nax2)
        HQ = CONJ(H)
        HABS = ABS(H)
        H2 = HABS^2
        GHQ = G*HQ  
       
        if keyword_set(psfout) then psfs[*,*,ngpsf] = psf
        ngpsf++

        nexp = nexp + thisfov
        GM = FFT(thisfov)
        GMHQ = GM*HQ
        masknumer = masknumer + GMHQ
        fimasknumer = fimasknumer + GMHQ

        numer = numer + GHQ
        denom = denom + H2
        nthis = nthis + 1
 
        finumer = finumer + GHQ
        fidenom = fidenom + H2
        ngood++
        
        n_skip = nthis mod nsub
        for j_sub = 0, nsub-1 do begin
          if (j_sub ne n_skip) then begin
            numersubreal[*,*,j_sub] = numersubreal[*,*,j_sub] + REAL_PART(GHQ)
            denomsubreal[*,*,j_sub] = denomsubreal[*,*,j_sub] + REAL_PART(H2)
            numersubim[*,*,j_sub] = numersubim[*,*,j_sub] + IMAGINARY(GHQ)
            denomsubim[*,*,j_sub] = denomsubim[*,*,j_sub] + IMAGINARY(H2)
            masknumersubreal[*,*,j_sub] = masknumersubreal[*,*,j_sub] + REAL_PART(GMHQ)
            masknumersubim[*,*,j_sub] = masknumersubim[*,*,j_sub] + IMAGINARY(GMHQ)
          endif
        endfor
        
        n_skip = ngood mod nsub
        for j_sub = 0, nsub-1 do begin
          if (j_sub ne n_skip) then begin
            finumersubreal[*,*,j_sub] = finumersubreal[*,*,j_sub] + REAL_PART(GHQ)
            fidenomsubreal[*,*,j_sub] = fidenomsubreal[*,*,j_sub] + REAL_PART(H2)
            finumersubim[*,*,j_sub] = finumersubim[*,*,j_sub] + IMAGINARY(GHQ)
            fidenomsubim[*,*,j_sub] = fidenomsubim[*,*,j_sub] + IMAGINARY(H2)
            fimasknumersubreal[*,*,j_sub] = fimasknumersubreal[*,*,j_sub] + REAL_PART(GMHQ)
            fimasknumersubim[*,*,j_sub] = fimasknumersubim[*,*,j_sub] + IMAGINARY(GMHQ)
          endif
        endfor

  
        hhh = REAL_PART(FFT(OTF*(GHQ/H2), /INVERSE))
        rcim = shift(hhh, center+1)

        if (debug gt 0) then begin
          mmm = REAL_PART(FFT(OTF*(GMHQ/H2), /INVERSE))
          rcmask = shift(mmm, center+1)
          writefits, tmpdir + 'rcmask.fits', rcmask
          writefits, tmpdir + 'rcim.fits', rcim
        endif
  
        print, "Frame number: " + string(i+1) + " ok"

       endif else begin
         print, "Frame number: " + string(i+1) + " no speckle cloud detected"; 
       endelse   ; if (count gt 1.0)


    ; output of control data
    ; ----------------------

    if ((nthis mod out_iter) eq 0 and nthis gt 0) then begin
     hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
     rcim = shift(hhh, center+1)
     mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
     rcmask = shift(mmm, center+1)
     support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
     rcim[support] = rcim[support]/rcmask[support]
     tmpexp = expmap
     if (nosup[0] gt -1) then begin
      tmpexp[nosup] = 0
      rcmask[nosup] = 0
      rcim[nosup] = 0
     endif
     writefits, tmpdir + 'rcim_tmp.fits', rcim
     writefits, tmpdir + 'rcmask_tmp.fits', rcmask
     writefits, tmpdir + 'expmap_tmp.fits', tmpexp
    endif

    if (debug gt 0) then STOP

    endfor                         ; end loop over frames in this cube

    endif else begin
      print, "Cube number: " + string(ic+1) + " no reference stars detected"
    endelse
  

  ; output of PSFs
  if keyword_set(psfout) then writefits, '../psfs/psfs' + strtrim(string(ic+1),2) + '.fits', psfs[*,*,0:ngpsf-1], /COMPRESS

  ; output of holography image for this cube
  ; ----------------------------------------------------
  hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
  rcim = shift(hhh, center+1)
  mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(expmap gt 0 and rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
  rcim[support] = rcim[support]/rcmask[support]
  if (nosup[0] gt -1) then begin
    expmap[nosup] = 0
    rcim[nosup] = 0
    rcmask[nosup] = 0
  endif
  writefits, tmpdir + list[ic] + '_support.fits', rcmask, /COMPRESS
  writefits, tmpdir + list[ic] + '_holo.fits', rcim, header, /COMPRESS
  writefits, tmpdir + list[ic] + '_expmap.fits', expmap, /COMPRESS
  printf, outfiles1, list[ic] + '_holo'
  printf, outfiles2, list[ic] + '_support'
  printf, outfiles3, list[ic] + '_expmap'

  ; output of sub-images and weight maps for this cube
  ; ---------------------------------------------------

  for is = 1, nsub do begin
   nhhh = complex(numersubreal[*,*,is-1],numersubim[*,*,is-1])
   dhhh = complex(denomsubreal[*,*,is-1],denomsubim[*,*,is-1])
   hhh = REAL_PART(FFT(OTF*(nhhh/dhhh), /INVERSE))
   rcim = shift(hhh, center+1)
   maskhhh = complex(masknumersubreal[*,*,is-1],masknumersubim[*,*,is-1])
   mmm = REAL_PART(FFT(OTF*(maskhhh/dhhh), /INVERSE))
   rcmask = shift(mmm, center+1)
   support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
   rcim[support] = rcim[support]/rcmask[support]
   if (nosup[0] gt -1) then begin
    rcim[nosup] = 0
    rcmask[nosup] = 0
   endif
   writefits, tmpdir + list[ic] + '_support' + '_s' + strtrim(string(is),2) +'.fits', (float(nthis)/nsub)*rcmask, /COMPRESS
   writefits, tmpdir + list[ic] + '_holo' + '_s' + strtrim(string(is),2) + '.fits', rcim, /COMPRESS
  endfor

 ; output of current total holography image
  ; ----------------------------------------------------
   hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
   rcim = shift(hhh, center+1)
   mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
   rcmask = shift(mmm, center+1)
   support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
   rcim[support] = rcim[support]/rcmask[support]
   tmpexp = nexp
   if (nosup[0] gt -1) then begin
    tmpexp[nosup] = 0
    rcim[nosup] = 0
    rcmask[nosup] = 0
   endif
   writefits, tmpdir + 'current_expmap.fits', tmpexp, /COMPRESS
   writefits, tmpdir + 'current_support.fits', nthis*rcmask, /COMPRESS
  writefits, tmpdir + 'current.fits', rcim, /COMPRESS

endfor ; end loop over input cubes


; final image
; --------------------
hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
rcim = shift(hhh, center+1)
mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
rcmask = shift(mmm, center+1)
support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
rcim[support] = rcim[support]/rcmask[support]
tmpexp = nexp
if (nosup[0] gt -1) then begin
 tmpexp[nosup] = 0
 rcim[nosup] = 0
 rcmask[nosup] = 0
endif
writefits, tmpdir + outnam+'_support.fits', nthis*rcmask, /COMPRESS
writefits, tmpdir + outnam + '_expmap.fits', tmpexp, /COMPRESS
writefits, tmpdir + outnam+'.fits', rcim, /COMPRESS
print, 'Total number of frames used: ' + strtrim(string(ngood),2)

; Output of image not re-convolved with OTF
if rawout gt 0 then begin
 hhh = REAL_PART(FFT(finumer/fidenom, /INVERSE))
 rcim = shift(hhh, center+1)
 mmm = REAL_PART(FFT(fimasknumer/fidenom, /INVERSE))
 rcmask = shift(mmm, center+1)
 support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
 rcim[support] = rcim[support]/rcmask[support]
 if (nosup[0] gt -1) then begin
  rcim[nosup] = 0.0
  rcmask[nosup] = 0
 endif
 writefits, tmpdir +  outnam+'_support_raw.fits', nthis*rcmask, /COMPRESS
 writefits, tmpdir + outnam+'_raw.fits', rcim, /COMPRESS
endif


; sub-images
; --------------------

for is = 1, nsub do begin
 nhhh = complex(finumersubreal[*,*,is-1],finumersubim[*,*,is-1])
 dhhh = complex(fidenomsubreal[*,*,is-1],fidenomsubim[*,*,is-1])
 hhh = REAL_PART(FFT(OTF*(nhhh/dhhh), /INVERSE))
 rcim = shift(hhh, center+1)
 maskhhh = complex(fimasknumersubreal[*,*,is-1],fimasknumersubim[*,*,is-1])
 mmm = REAL_PART(FFT(OTF*(maskhhh/dhhh), /INVERSE))
 rcmask = shift(mmm, center+1)
 support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
 rcim[support] = rcim[support]/rcmask[support]
 if (nosup[0] gt -1) then begin
  rcim[nosup] = 0
  rcmask[nosup] = 0
 endif
 writefits, tmpdir + outnam + '_support_s' + strtrim(string(is),2) + '.fits',  (float(nthis)/nsub)*rcmask, /COMPRESS
 writefits, tmpdir + outnam + '_s' + strtrim(string(is),2) + '.fits', rcim, /COMPRESS
endfor

print, 'Number of frames used in sub-images: ' + strtrim(string(float(ngood/nsub)),2)

; Output of images not re-convolved with OTF
;if rawout gt 0 then begin
; for is = 1, nsub do begin
;  nhhh = complex(finumersubreal[*,*,is-1],finumersubim[*,*,is-1])
;  dhhh = complex(fidenomsubreal[*,*,is-1],fidenomsubim[*,*,is-1])
;  hhh = REAL_PART(FFT(nhhh/dhhh, /INVERSE))
;  rcim = shift(hhh, center+1)
;  maskhhh = complex(fimasknumersubreal[*,*,is-1],fimasknumersubim[*,*,is-1])
;  mmm = REAL_PART(FFT(maskhhh/dhhh, /INVERSE))
;  rcmask = shift(mmm, center+1)
;  support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
;  rcim[support] = rcim[support]/rcmask[support]
;  if (nosup[0] gt -1) then begin
;   rcim[nosup] = 0
;   rcmask[nosup] = 0
;  endif
;  writefits, outnam + '_s_raw' + strtrim(string(is),2) + '.fits', rcim, /COMPRESS
; endfor
;endif

free_lun, outfiles1, outfiles2, outfiles3

END
 
