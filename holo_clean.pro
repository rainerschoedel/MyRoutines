PRO CALEVAL, X, P, YMOD

  YMOD = P[0] + P[1]*X
  YMOD = image_shift(YMOD,P[2],P[3])

END


;---------------------

PRO HOLO_CLEAN, indir, inlist, outnam,  maskrad, normrad, nsigma, rebfac, refsources, starlist, maskdir=maskdir, masklist=masklist, CLEANDIR = cleandir, DEBUG = debug, iter=iter, AIRY=airy, OUT_ITER=out_iter, PSFOUT = psfout, UNWEIGHTED = unweighted, SUBPIX=subpix, tmpdir=tmpdir, BOXHW=boxhw, NSUB=nsub, MINSUPP=minsupp, MAXSUPP = maxsupp, CMINSUPP=cminsupp, CMAXSUPP = cmaxsupp,   REBITER=rebiter, RAWOUT = rawout, N_MASK_SECONDARY=n_mask_secondary, PSFNOISE = psfnoise, PSF_FRAC = psf_frac, CORRECT_SKY = correct_sky, SMOOTHMASK = smoothmask, N_REF_MAX = n_ref_max, PR = pr, CIRC_BORDER=circ_border, CleanScale = CleanScale

; ---------------------------------------------------
; As HOLO_FOURIER, but subtracts the noise threshold
; ---------------------------------------------------

withmask = KEYWORD_SET(masklist)
if (not KEYWORD_SET(debug)) then debug = 0
if (not KEYWORD_SET(maxsupp)) then maxsupp = 1.1
if (not KEYWORD_SET(minsupp)) then minsupp = 0.01
if (not KEYWORD_SET(cmaxsupp)) then cmaxsupp = 1.1
if (not KEYWORD_SET(cminsupp)) then cminsupp = 0.9
if (not KEYWORD_SET(n_mask_secondary)) then n_mask_secondary = 0
if (not KEYWORD_SET(rawout)) then rawout= 0
if (not KEYWORD_SET(psfnoise)) then psfnoise = 0
if (not KEYWORD_SET(scale_filter)) then scale_filter = 1.e-12
if (not KEYWORD_SET(rebiter)) then rebiter = 0
if (not KEYWORD_SET(psf_frac)) then psf_frac = 0.9
if (not KEYWORD_SET(correct_sky)) then correct_sky = 0
if (not KEYWORD_SET(smoothmask)) then smoothmask = 0
if (not KEYWORD_SET(pr)) then pr = 1
if (not KEYWORD_SET(circ_border)) then circ_border = 0
if (not KEYWORD_SET(CleanScale)) then CleanScale = 4


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
dummy = fltarr(nax1,nax2)
dummy[*,*] = 1
bgring0 = circ_mask(dummy,cen1,cen2,maskrad, INNER=1)
bgring0 = circ_mask(bgring0,cen1,cen2,2*maskrad)
bgind0 = where(bgring0 gt 0, n_bgind0)

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

fimask = fltarr(nax1,nax2)
finumer = fltarr(nax1,nax2)
fidenom = fltarr(nax1,nax2)
finumersubreal = fltarr(nax1,nax2,nsub)
fidenomsubreal = fltarr(nax1,nax2,nsub)
finumersubim = fltarr(nax1,nax2,nsub)
fidenomsubim = fltarr(nax1,nax2,nsub)
ngood = 0L

if (withmask) then begin
 fimasknumer = fltarr(nax1,nax2)
 fimasknumersubreal = fltarr(nax1,nax2,nsub)
 fimasknumersubim = fltarr(nax1,nax2,nsub)
endif

; reference sources
x_psf = refsources[0,*]
y_psf = refsources[1,*]
; rebin positions if neessary
if (rebiter lt 1) then begin
 x_psf = rebfac * x_psf
 y_psf = rebfac * y_psf
endif
nref = n_elements(x_psf)
xint = round(x_psf)
yint = round(y_psf)

; open files for output lists
openw, outfiles1, tmpdir + 'holo.txt', /get_lun
openw, outfiles2, tmpdir + 'support.txt', /get_lun

; delta function
delta_width = 7*pr
delta = fltarr(delta_width,delta_width)
delta_cen = delta_width/2
delta[delta_cen,delta_cen] = 1.0

for ic = 0, ncubes-1 do begin  ; start loop over all input cubes

  ; initialize variables for holography
  ;---------------------------------------

  expmap = replicate(0.0,nax1,nax2)
  mask = replicate(0.0,nax1,nax2)
;  mask_sub = replicate(0.0,nax1,nax2,nsub)
  masknumer = replicate(0.0,nax1,nax2)
  numer = replicate(0.0,nax1,nax2)
  denom = replicate(0.0,nax1,nax2)
  numersubreal = fltarr(nax1,nax2,nsub)
  denomsubreal = fltarr(nax1,nax2,nsub)
  numersubim = fltarr(nax1,nax2,nsub)
  denomsubim = fltarr(nax1,nax2,nsub)
  if (withmask) then begin
   masknumersubreal = fltarr(nax1,nax2,nsub)
   masknumersubim = fltarr(nax1,nax2,nsub)
  endif
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

;  rcims = fltarr(nax1,nax2,nim)

  if keyword_set(psfout) then psfs = fltarr(nax1,nax2,nim)
  ngpsf = 0L


  for i = 0L, nim-1 do begin  ; start loop over frames in this cube

    ; load ith image 
    im = cube[*,*,i]

   if (withmask) then begin
    thisfov = mcube[*,*,i]  
    im = im * thisfov ; in case im has non-zero values outside of FOV
   endif else thisfov = replicate(1,nax1,nax2)

   ; Rebin image and mask, smooth mask edges (even if rebfac = 1)
   im = CREBIN(im,nax1,nax2,/TOTAL)
   if (withmask) then begin
    thisfov = CREBIN(thisfov,nax1,nax2) ; Do not normalize rebinned FOV!
    edgpix = where(thisfov ne 1)
    if edgpix[0] gt -1 then thisfov[edgpix] = 0
    im = im * thisfov
    expmap = expmap + thisfov
   endif
 
  ; PSF extraction
  ; ################

    ; select valid reference sources
    valid_ref = replicate(1,nref)
    sub_arrays, im, xint, yint, boxsize, stack, masks
    if withmask then begin
     sub_arrays, thisfov, xint, yint, boxsize, psfmasks, masks
     for iref = 0, nref-1 do begin
       masks[*,*,iref] = masks[*,*,iref]*psfmasks[*,*,iref]
       dummy = circ_mask(masks[*,*,iref],boxhw,boxhw,maskrad)
       if total(dummy) lt round(psf_frac*npix_ref) then begin
         valid_ref[iref] = 0
       endif
     endfor
    endif
    acceptref = where(valid_ref gt 0, nref_accept,complement=reject)
    xint_accept = xint[acceptref]
    yint_accept = yint[acceptref]
    x_psf_accept = x_psf[acceptref]
    y_psf_accept = y_psf[acceptref]
    xint_reject = xint[reject]
    yint_reject = yint[reject]

    tmpmask = thisfov
    tmpmask = CREBIN(tmpmask,pr*nax1,pr*nax2) ; Do not normalize rebinned FOV!
    edgpix = where(tmpmask ne 1)
    if edgpix[0] gt -1 then tmpmask[edgpix] = 0
    mask_reject = tmpmask
    for j_reject = 0, n_elements(reject)-1 do begin
     mask_reject = mask_reject * circ_mask(tmpmask,pr*xint_reject[j_reject],pr*yint_reject[j_reject],pr*maskrad/2,/INNER)
    endfor

    
    delta_im = image_model(xx*pr,yy*pr,ff,nax1*pr,nax2*pr,delta,REFERENCE_PIX=[delta_cen,delta_cen])
    delta_im = delta_im*mask_reject
;    delta_mask = image_model(xx*pr,yy*pr,1,nax1*pr,nax2*pr,delta,REFERENCE_PIX=[delta_cen,delta_cen])
;    delta_mask = delta_mask*mask_reject

    psfim = CREBIN(im,pr*nax1,pr*nax2,/TOTAL)
    
    G = FFT(psfim*mask_reject)
    H = FFT(delta_im)*(nax1*nax2)
    GH = G/H
    GH_REAL = REAL_PART(GH)
    GH_REAL = sigma_filter(GH_REAL, 3*pr, N_SIGMA=5., /ITERATE, /ALL_PIXELS)
    GH_IM   = IMAGINARY(GH)
    GH_IM = sigma_filter(GH_IM, 3*pr, N_SIGMA=5., /ITERATE, /ALL_PIXELS)
    GH = COMPLEX(GH_REAL, GH_IM)


    hhh = REAL_PART(FFT(GH,/INVERSE))
    hhh = CREBIN(hhh,nax1,nax2,/TOTAL)
    psf = shift(hhh, center+1)
    
     ; mask PSF
     ; -----------------------
     if (debug gt 0) then begin
      writefits, tmpdir + 'im.fits', im
      if withmask then writefits, tmpdir + 'thisfov.fits', thisfov
      writefits, tmpdir + 'bgring_0.fits', psf*bgring
;      writefits, tmpdir + 'psfmod.fits', psf_mask
      writefits, tmpdir + 'psfmask.fits', mask_reject
      writefits, tmpdir + 'psf_0raw.fits', psf
     endif
;      print, 'Background: ' + string(median(psf[bgind]))
;      use astrolib robus statistics script
;      to estimate mean and stddev
;      this helps to avoid rare problems with individual bright pixels
      RESISTANT_Mean,psf[bgind0],3.0,mean_psf,sigma_psf,Num_Rej
      noise = sigma_psf * sqrt(n_bgind0 - Num_Rej)
      psf = psf - mean_psf- nsigma[0] * noise
      suppress = where(psf lt 0, count)        
;      suppress = where(psf lt nsigma[0] * noise, count)        
      if (count gt 0) then psf[suppress] = 0

      sub_arrays, psf, cen1, cen2, boxsize, slice
      psf = slice
      MASK_PSF, psf, maskrad, PSF_MASKED=psf_masked
      psf = psf_masked
      psf = psf/total(psf)
      psf_fwhm = round(fwhm(psf))
      if (debug gt 0) then begin
        writefits, tmpdir + 'deltaim.fits', delta_im
        writefits, tmpdir + 'psf_0.fits', psf
      endif
    
 

  ; improval of PSF with known sources
  ; ####################################

      for it = 0, iter-1 do begin

       ; clean surroundings of reference stars from secondary sources
       ; ------------------------------------------------------------
      stack = fltarr(boxsize,boxsize,nref_accept)
      psfmasks = fltarr(boxsize,boxsize,nref_accept)
      skyarr = fltarr(nref_accept) ; to store sky values

                                ; use psfmasks if a reference star may not be contained in all
      ; cubes because of dithering
      ; ---------------------------

      for iref = 0, nref_accept-1 do begin
       sub_arrays, im, xint_accept[iref], yint_accept[iref], boxsize, slice, masks
       if withmask then begin
         sub_arrays, thisfov, xint_accept[iref], yint_accept[iref], boxsize, maskslice, masks
         psfmasks[*,*,iref] = maskslice * masks
        endif else psfmasks[*,*,iref] = masks
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
         model = image_model(x_near,y_near,f_near,boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
         P = [0.,1.,0.,0.]
         W = model * psfmasks[*,*,iref]
         res = mpcurvefit(model,slice,W,P,sigma,FUNCTION_NAME='CALEVAL',/NODERIVATIVE,/QUIET)
         compare_lists, x_near, y_near, x_psf_accept[iref]-(xint_accept[iref]-boxhw), y_psf_accept[iref]-(yint_accept[iref]-boxhw), x1c, y1c, x2c ,y2c, SUB1=other_stars, MAX_DISTANCE=dmax 
         if (other_stars[0] gt -1) then begin
          secondaries = image_model(x_near[other_stars],y_near[other_stars],f_near[other_stars],boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
;        if (debug gt 0) then begin
;         print, P
;         print, iref
;         writefits, tmpdir + 'model.fits', model
;         writefits, tmpdir + 'rawslice.fits', slice
;         writefits, tmpdir + 'secondaries.fits', image_shift(secondaries*P[1],P[2],P[3])
;         writefits, tmpdir + 'diff.fits', slice - image_shift(secondaries*P[1],P[2],P[3])
;         STOP
;        endif
         slice = slice - image_shift(P[0]+secondaries*P[1],P[2],P[3])
         skyarr[iref] = P[0]
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
           psf_mask = circ_mask(psf,xymax[0],xymax[1],round(psf_FWHM))
           x_mask = round((x_near[other_stars[0:n_mask-1]]))
           y_mask = round((y_near[other_stars[0:n_mask-1]]))
           f_mask = round((f_near[other_stars[0:n_mask-1]]))  ; just needed for formal reasons
           mask_secondaries = image_model(x_mask,y_mask,f_mask,boxsize,boxsize,psf_mask,REFERENCE_PIX=[boxhw,boxhw])
           mask_ind = where(mask_secondaries gt 0, complement=nomask)
           maskslice = psfmasks[*,*,iref]
           if (mask_ind[0] gt -1) then maskslice[mask_ind] = 0
            psfmasks[*,*,iref] = maskslice
           endif      ;  if (n_mask_secondary gt 0)
          endif       ;  (other_stars[0] gt -1)
        endif        ;  if (nearby[0] gt -1)
        stack[*,*,iref] = slice
       endfor
       sky = mean(skyarr)
      
       if (debug gt 0) then writefits,  tmpdir + 'rawstack_1.fits', stack
       for iref = 0, nref_accept-1 do begin
         subim = stack[*,*,iref]
         submask = psfmasks[*,*,iref]
;        sky_pixels = subim[bgind]
;        sky_mask = submask[bgind]
;        sky_good = where(sky_mask gt 0)
;        bg_level = median(sky_pixels[sky_good])
;        sky = sky + bg_level
;        print, bg_level
         ; constant background should already be subtracted (P[0])
;        normim = subim - bg_level
;        print, 'sky for reference star ' + strn(iref+1) + ':  ' + strn(bg_level)
        ; the next IF statement avoids a problem
        ; if the entire subim = 0 (can occur when a reference
        ; star lies outside the current FOV
        normim = subim
        if min(normim lt 0) then normim[where(normim lt 0)] = 0
        normim = circ_mask(normim, boxhw, boxhw, normrad)
        if (unweighted eq 1 and total(normim) gt 0) then begin
;         subim = subim/total(normim)
          subim = subim/max(normim)
        endif
        if (subpix gt 0) then begin
         subim = image_shift(subim,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
         submask = image_shift(submask,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
         zeroind = where(submask lt 0.999, complement = ones)
         if (zeroind[0] gt -1) then begin
          submask[zeroind] = 0
          submask[ones] = 1
         endif  
        endif
        stack[*,*,iref] = subim
        psfmasks[*,*,iref] = submask
       endfor
       psf = stack_median(stack, MASK=psfmasks)
       psf_sigma = stack_error(stack, MASK=psfmasks)/sqrt(nref_accept)

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
       psf = psf - mean_bg
;      Noise suppression with background ring if there are
;      less than 3 PSF reference stars or if psfnoise  < 1
       if (nref_accept lt 3) or (psfnoise lt 1) then begin
          psf = psf - nsigma[it+1] * noise
          suppress = where(psf  lt 0, complement=accept,count)
       endif else begin
;      Noise suppression with PSF noise
         psf = psf - nsigma[it+1] * psf_sigma
         suppress = where(psf lt 0, complement=accept,count)    
       endelse  
       if (count gt 0) then psf[suppress] = 0
       ;psfcen = centroid(circ_mask(psf,boxhw,boxhw,maskrad))
       psf = circ_mask(psf,boxhw,boxhw,maskrad,BORDER=circ_border)
       psfnorm = total(psf)
       psf = psf/psfnorm
       if (debug gt 0) then begin
        writefits, tmpdir + 'stack.fits', stack*psfmasks
        writefits, tmpdir + 'stackmasks.fits', psfmasks
        writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + '.fits', psf
     endif
     endfor ; end loop over iter

    ; clean image from artifacts
     ; ##########################

     model =  image_model(xx,yy,ff,nax1,nax2,psf,REFERENCE_PIX=[boxhw,boxhw])
     P = [0.,1.,0.,0.]
     W = thisfov
     res = mpcurvefit(model,im,W,P,sigma,FUNCTION_NAME='CALEVAL',/NODERIVATIVE,/QUIET)
     model = image_shift(P[0]+model*P[1],P[2],P[3])
     resid = im - model
     clean = del_pattern(resid,pattern=pattern,NbrScale=CleanScale,NSigma=5)
     if (debug gt 0) then begin
       writefits, tmpdir + 'resid.fits', resid
       writefits, tmpdir + 'resid_clean.fits', clean
       writefits, tmpdir + 'clean.fits', clean + model
       writefits, tmpdir + 'pattern.fits', pattern
     endif
     im = clean + model
     cube[*,*,i] = CREBIN(im,nax1/rebfac,nax2/rebfac,/TOTAL)

     ; adapt PSF to size of input frame
     ; shift it to the center (only integer shift!)
     temp = fltarr(nax1,nax2)
     temp[0:boxsize-1,0:boxsize-1] = psf
     psf = shiftnw(temp,cen1-boxhw,cen2-boxhw)

 

; holography algorithm
  ; #####################

      ind = where(psf gt 0.0, count)

      if (count gt 0)  then begin
         
        if (correct_sky gt 0) then begin
         print, 'sky from PSF: ' + strn(sky)
         im = thisfov*(im - sky)  ; subtract sky
        endif
 

      ; smooth mask at edges
       if (withmask) then begin
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
        im = im * thisfov  ; for correct weighting of mask edges
 ;     writefits, tmpdir + 'thisfov.fits', thisfov
 ;     writefits, tmpdir + 'im.fits', im
 ;     STOP
       endif

        G = FFT(im)
        H = FFT(psf)*(nax1*nax2)
        HQ = CONJ(H)
        HABS = ABS(H)
        H2 = HABS^2
        GHQ = G*HQ
  
       
        if keyword_set(psfout) then psfs[*,*,ngpsf] = psf
        ngpsf++

        if (withmask) then begin
         mask = mask + thisfov
         fimask = fimask + thisfov
         GM = FFT(thisfov)
         GMHQ = GM*HQ
         masknumer = masknumer + GMHQ
         fimasknumer = fimasknumer + GMHQ
        endif

        numer = numer + GHQ
        denom = denom + H2
        nthis = nthis + 1
 
        finumer = finumer + GHQ
        fidenom = fidenom + H2
        ngood++
        
        n123 = nthis mod nsub
        numersubreal[*,*,n123] = numersubreal[*,*,n123] + REAL_PART(GHQ)
        denomsubreal[*,*,n123] = denomsubreal[*,*,n123] + REAL_PART(H2)
        numersubim[*,*,n123] = numersubim[*,*,n123] + IMAGINARY(GHQ)
        denomsubim[*,*,n123] = denomsubim[*,*,n123] + IMAGINARY(H2)
        if (withmask) then begin
;         mask_sub[*,*,n123] = mask_sub[*,*,n123] + thisfov
         masknumersubreal[*,*,n123] = masknumersubreal[*,*,n123] + REAL_PART(GMHQ)
         masknumersubim[*,*,n123] = masknumersubim[*,*,n123] + IMAGINARY(GMHQ)
        endif
        
        n123 = ngood mod nsub
        finumersubreal[*,*,n123] = finumersubreal[*,*,n123] + REAL_PART(GHQ)
        fidenomsubreal[*,*,n123] = fidenomsubreal[*,*,n123] + REAL_PART(H2)
        finumersubim[*,*,n123] = finumersubim[*,*,n123] + IMAGINARY(GHQ)
        fidenomsubim[*,*,n123] = fidenomsubim[*,*,n123] + IMAGINARY(H2)
        if (withmask) then begin
         fimasknumersubreal[*,*,n123] = fimasknumersubreal[*,*,n123] + REAL_PART(GMHQ)
         fimasknumersubim[*,*,n123] = fimasknumersubim[*,*,n123] + IMAGINARY(GMHQ)
        endif

  
        hhh = REAL_PART(FFT(OTF*(GHQ/H2), /INVERSE))
        rcim = shift(hhh, center+1)
;        rcims[*,*,i] = rcim

        if (debug gt 0) then begin
           if withmask then begin
              mmm = REAL_PART(FFT(OTF*(GMHQ/H2), /INVERSE))
              rcmask = shift(mmm, center+1)
              writefits, tmpdir + 'rcmask.fits', rcmask
            endif
            writefits, tmpdir + 'rcim.fits', rcim
        endif
  
        print, "Frame number: " + string(i+1) + " ok"

       endif else begin
         print, "Frame number: " + string(i+1) + " no speckle cloud detected"
       endelse   ; if (count gt 1.0)

       if (debug gt 0) then STOP



    ; output of control data
    ; ----------------------

    if ((nthis mod out_iter) eq 0 and nthis gt 0) then begin
     hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
     rcim = shift(hhh, center+1)
     if (withmask) then begin
      mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
      rcmask = shift(mmm, center+1)
      writefits, tmpdir + 'rcmask.fits', rcmask
      support = where(expmap gt cminsupp*nthis and rcmask lt cmaxsupp*nthis, complement=nosup)
      rcim[support] = rcim[support]/rcmask[support]
      if (nosup[0] gt -1) then begin
       rcim[nosup] = 0
      endif
     endif
     writefits, tmpdir + 'rcim_tmp.fits', rcim
    endif
   endfor  ; end loop over frames in this cube
  
   ; Write clean images to disk
   ; ###########################
   writefits, cleandir + list[ic], cube, header


  ; output of PSFs
  if keyword_set(psfout) then writefits, '../psfs/psfs' + strtrim(string(ic+1),2) + '.fits', psfs[*,*,0:ngpsf-1]

  ; output of holography image for this cube
  ; ----------------------------------------------------
  hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
  rcim = shift(hhh, center+1)
  if (withmask) then begin
   mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
   rcmask = shift(mmm, center+1)
   support = where(expmap gt cminsupp*nthis and rcmask lt cmaxsupp*nthis, complement=nosup)
   if (support[0] gt -1) then begin
    rcim[support] = rcim[support]/rcmask[support]
    if (nosup[0] gt -1) then begin
     rcim[nosup] = 0
     rcmask[nosup] = 0
    endif
   endif
   writefits, tmpdir + list[ic] + '_mask.fits', mask
   writefits, tmpdir + list[ic] + '_support.fits', nthis*rcmask
  endif
  writefits, tmpdir + list[ic] + '_holo.fits', rcim, header
  writefits, tmpdir + list[ic] + '_expmap.fits', expmap
  printf, outfiles1, list[ic] + '_holo'
  printf, outfiles2, list[ic] + '_support'

  ; output of sub-images and weight maps for this cube
  ; ---------------------------------------------------

  for is = 1, nsub do begin
   nhhh = complex(numersubreal[*,*,is-1],numersubim[*,*,is-1])
   dhhh = complex(denomsubreal[*,*,is-1],denomsubim[*,*,is-1])
   hhh = REAL_PART(FFT(OTF*(nhhh/dhhh), /INVERSE))
   rcim = shift(hhh, center+1)
   if (withmask) then begin
    maskhhh = complex(masknumersubreal[*,*,is-1],masknumersubim[*,*,is-1])
    mmm = REAL_PART(FFT(OTF*(maskhhh/dhhh), /INVERSE))
    rcmask = shift(mmm, center+1)
    support = where(rcmask gt cminsupp and rcmask lt cmaxsupp, complement=nosup)
    rcim[support] = rcim[support]/rcmask[support]
    if (nosup[0] gt -1) then begin
     rcim[nosup] = 0
     rcmask[nosup] = 0
    endif
    writefits, tmpdir + list[ic] + '_support' + '_s' + strtrim(string(is),2) +'.fits', (float(nthis)/nsub)*rcmask
   endif
   writefits, tmpdir + list[ic] + '_holo' + '_s' + strtrim(string(is),2) + '.fits', rcim
  endfor

 ; output of current total holography image
  ; ----------------------------------------------------
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
   writefits, tmpdir + 'current_mask.fits', fimask
   writefits, tmpdir + 'current_support.fits', ngood*rcmask
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
 writefits, outnam+'_support.fits', ngood*rcmask
 writefits, outnam + '_mask.fits', fimask
endif
writefits, outnam+'.fits', rcim
print, 'Total number of frames used: ' + strtrim(string(ngood),2)

; Output of image not re-convolved with OTF
if rawout gt 0 then begin
 hhh = REAL_PART(FFT(finumer/fidenom, /INVERSE))
 rcim = shift(hhh, center+1)
 if (withmask) then begin
  mmm = REAL_PART(FFT(fimasknumer/fidenom, /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
  rcim[support] = rcim[support]/rcmask[support]
  if (nosup[0] gt -1) then begin
   rcim[nosup] = 0.0
   rcmask[nosup] = 0
  endif
  writefits, outnam+'_support_raw.fits', ngood*rcmask
  writefits, outnam + '_mask_raw.fits', fimask
 endif
 writefits, outnam+'_raw.fits', rcim
endif


; sub-images
; --------------------

for is = 1, nsub do begin
 nhhh = complex(finumersubreal[*,*,is-1],finumersubim[*,*,is-1])
 dhhh = complex(fidenomsubreal[*,*,is-1],fidenomsubim[*,*,is-1])
 hhh = REAL_PART(FFT(OTF*(nhhh/dhhh), /INVERSE))
 rcim = shift(hhh, center+1)
 if (withmask) then begin
  maskhhh = complex(fimasknumersubreal[*,*,is-1],fimasknumersubim[*,*,is-1])
  mmm = REAL_PART(FFT(OTF*(maskhhh/dhhh), /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
  rcim[support] = rcim[support]/rcmask[support]
  if (nosup[0] gt -1) then begin
   rcim[nosup] = 0
   rcmask[nosup] = 0
  endif
 endif
 writefits, outnam + '_s' + strtrim(string(is),2) + '.fits', rcim
endfor

print, 'Number of frames used in sub-images: ' + strtrim(string(float(ngood/nsub)),2)

; Output of images not re-convolved with OTF
if rawout gt 0 then begin
 for is = 1, nsub do begin
  nhhh = complex(finumersubreal[*,*,is-1],finumersubim[*,*,is-1])
  dhhh = complex(fidenomsubreal[*,*,is-1],fidenomsubim[*,*,is-1])
  hhh = REAL_PART(FFT(nhhh/dhhh, /INVERSE))
  rcim = shift(hhh, center+1)
  if (withmask) then begin
   maskhhh = complex(fimasknumersubreal[*,*,is-1],fimasknumersubim[*,*,is-1])
   mmm = REAL_PART(FFT(maskhhh/dhhh, /INVERSE))
   rcmask = shift(mmm, center+1)
   support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
   rcim[support] = rcim[support]/rcmask[support]
   if (nosup[0] gt -1) then begin
    rcim[nosup] = 0
    rcmask[nosup] = 0
   endif
  endif
  writefits, outnam + '_s_raw' + strtrim(string(is),2) + '.fits', rcim
 endfor
endif

free_lun, outfiles1, outfiles2

END
 
