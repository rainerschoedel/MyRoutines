;---------------------
; MODIFICATION HISTORY
;
; 24 July 2024, Rainer Schoedel
; USE_CENTROID can be set to sub-pixel shift the reference stars with the 
; StarFinder CENTROID function. This is necessary if the positions
; of the reference stars have not yet been determined with high accuracy.
; --------------------

PRO SSA_JACK, indir, inlist, out_path, outnam,  mrad, normrad, nsigma, refsources, rebfac, maskdir=maskdir, masklist=masklist, DEBUG = debug, UNWEIGHTED = unweighted,  tmpdir=tmpdir,  NJACK=njack, PSF_FRAC = psf_frac, CORRECT_SKY = correct_sky,  N_REF_MAX = n_ref_max, CIRC_BORDER=circ_border, SATLEVEL=satlevel, PSF_BORDER=psf_border, ESTIM_BG = estim_bg, SAT_MASK = sat_mask, contrast_thresh = contrast_thresh, USE_CENTROID=use_centroid, CEN_SHIFT=cen_shift, ROBUST=robust, BOX_HWIDTH = box_hwidth

; ---------------------------------------------------
; NOTE: cubes with less than njack frames will be discarded
; ---------------------------------------------------

withmask = KEYWORD_SET(masklist)
if (not KEYWORD_SET(debug)) then debug = 0
if (not KEYWORD_SET(psf_frac)) then psf_frac = 0.9
if (not KEYWORD_SET(correct_sky)) then correct_sky = 0
if (not KEYWORD_SET(circ_border)) then circ_border = 0
if (not KEYWORD_SET(satlevel)) then satlevel = 1.0e9
if (not KEYWORD_SET(psf_border)) then psf_border = 0.0
if (not KEYWORD_SET(estim_bg)) then estim_bg = 0
if (not KEYWORD_SET(contrast_thresh)) then contrast_tresh = 0
 if not(keyword_set(use_centroid)) then use_centroid = 0


; Identical to holo7sub.pro except that it 
; subtracts the noise threshold from the final PSF
; -------------------------------------------------

; read list of input cubes
readcol, inlist, list, FORMAT='A'
ncubes = n_elements(list)
print, 'input file list: ' 
print, list
print

; read list of  FOV masks
if (withmask) then begin
  readcol, masklist, mlist, FORMAT='A'
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
maskrad = rebfac * mrad
nax1 = rebfac*nax1
nax2 = rebfac*nax2
normrad = rebfac*normrad
boxhw = rebfac*box_hwidth ; change name of variable
; pixel indixes for subarray extraction
boxsize = long(2*boxhw + 1)
cen = boxsize/2
;cen1 = nax1/2
;cen2 = nax2/2
;center = [cen1,cen2]

; dummy needed  ringmask for background estimation
; npix_ref if the minimum required support region for
; reference stars (in case of variable masking/dithering)
dummy = fltarr(boxsize,boxsize)
dummy[*,*] = 1
bgring = circ_mask(dummy,boxhw,boxhw,maskrad, INNER=1)
bgind = where(bgring gt 0, n_bgind)
dummy[*,*] = 1
dummy = circ_mask(dummy,boxhw,boxhw,maskrad)
psf_support = where(dummy gt 0, npix_ref)

; ######################################
; reference sources
x_psf = refsources[0,*] * rebfac
y_psf = refsources[1,*] * rebfac
f_psf = refsources[2,*]/rebfac^2
nref = n_elements(x_psf)
xint = round(x_psf)
yint = round(y_psf)




; initialize variables for SSA
;---------------------------------------
ssa_all = fltarr(nax1,nax2)
ssa_sigma_all = fltarr(nax1,nax2)
wt_all = fltarr(nax1,nax2)
ssa_jack = fltarr(nax1,nax2,njack)
wt_jack = fltarr(nax1,nax2,njack)
nt = 0L  ; total number of frames used
ngood = 0L
nthis = 0L

; Now start lop over all input cubes
; --------------------------------

for ic = 0, ncubes-1 do begin


   if (withmask) then begin
    mcube = readfits(maskdir + mlist[ic])
   endif else mcube = replicate(1,nax1,nax2) ; define thisfov for exposure map
   fov_tot = total(long(mcube),3)
   bad = where(fov_tot lt max(fov_tot),complement=good)
   fov_tot[good] = 1
   fov_tot[bad] = 0
   if (debug gt 0) then begin
       writefits, tmpdir + 'fov_tot.fits', fov_tot
   endif
   fov_tot = CONGRID(fov_tot,nax1,nax2,/CENTER,CUBIC = -0.5)
   bad = where(fov_tot lt 1,complement=good)
   fov_tot[good] = 1
   fov_tot[bad] = 0

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
  sub_arrays, fov_tot, xint, yint, boxsize, psffovs, masks
  for iref = 0, nref-1 do begin
    tmpmask = masks[*,*,iref]*psffovs[*,*,iref]
    dummy = circ_mask(tmpmask,boxhw,boxhw,maskrad)
    if total(dummy) lt round(psf_frac*npix_ref) then begin
      valid_ref[iref] = 0
    endif
  endfor
  acceptref = where(valid_ref gt 0, nref_accept,complement=reject)

  if (nref_accept gt 0) and (nim ge njack) then begin
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

 
    if keyword_set(psfout) then begin
       psfs = fltarr(nax1,nax2,nim)
       psf_sigmas = fltarr(nax1,nax2,nim)      
    endif
    ngpsf = 0L

  ; initialize images
  ssa = fltarr(nax1,nax2)
  wt = fltarr(nax1,nax2)
  ssa_sigma = fltarr(nax1,nax2)

  outcube = fltarr(nax1,nax2,nim)  
  outmasks = fltarr(nax1,nax2,nim)  

  ; start loop over frames in this cube
  ; ----------------------------------
  for i = 0L, nim-1 do begin

     ; load ith image 
    im = cube[*,*,i]
    mask = mcube[*,*,i] 

    ; optionally rebin image and mask
    if (rebfac gt 1) then begin
      im = CONGRID(im,nax1,nax2,/CENTER,CUBIC = -0.5)
      mask = CONGRID(mask,nax1,nax2,/CENTER,CUBIC = -0.5)
      reject = where(mask lt 0.99, complement=accept)
      if (reject[0] gt -1) then begin
         mask[reject] = 0
         mask[accept] = 1
         im = im * mask
      endif
    endif
    thisfov = mask 

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
      if (use_centroid) then begin
        subim = CENTROIDER(subim,XSHIFT=xs,YSHIFT=ys)
        submask = image_shift(submask,xs,ys)
      endif else begin
        subim = image_shift(subim,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
      endelse
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
     psf = psf - mean_psf- nsigma * noise
     suppress = where(psf lt 0, count)        
     if (count gt 0) then psf[suppress] = 0
     psfcen = centroid(circ_mask(psf,boxhw,boxhw,maskrad))
     psf = circ_mask(psf,psfcen[0],psfcen[1],maskrad)
;     psf = circ_mask(psf,boxhw,boxhw,maskrad)
     psf = psf/total(psf)
     ; compute PSF FWHM for later use
     psf_fwhm = round(fwhm(psf))
     psf_sigma[*,*] = noise
     contrast = max(psf)/max(psf_sigma)
     print, 'PSF contrast: '+ strn(contrast)

     if (debug gt 0) then begin
       writefits, tmpdir + 'im.fits', im
       writefits, tmpdir + 'thisfov.fits', thisfov
       writefits, tmpdir + 'stack0.fits', stack
       writefits, tmpdir + 'bgring_0.fits', psf*bgring
       writefits, tmpdir + 'psf_0.fits', psf
    endif


    psfmax = max(psf, max_subscript)
    print, string(i+1) + ' ' + string(psfmax)
    mind = array_indices(psf,max_subscript)
    if (debug gt 0) then begin
      writefits, tmpdir + 'psf_preshift.fits', psf
    endif    

; shift each image and mask to position defined by PSF maximum or centroid
; ########################################################################

    neg = where(psf lt 0, complement=positive)
    if (neg[0] gt -1) then psf[neg] = 0
    if (cen_shift gt 0) then begin
     ; Do a modified centroid shift,
     ; i.e. shift to center of non-zero pixels
     cenpsf = psf
     cenpsf[positive] = 1.0
     centr = centroid(cenpsf)
     xs = round(cen - centr[0])
     ys = round(cen - centr[1])
    endif else begin
     xs = round(cen - mind[0])
     ys = round(cen - mind[1])
    endelse
    outcube[*,*,nthis] = shiftnw(im,xs,ys)
    outmasks[*,*,nthis] = shiftnw(mask,xs,ys)
    if (debug gt 0) then begin
       print, 'shift: ' + string(xs) + ', ' + string(ys)
       writefits, tmpdir + 'psf_postshift.fits', shiftnw(psf,xs,ys)
       writefits, tmpdir + 'mask.fits', outmasks[*,*,nthis]
       writefits, tmpdir + 'im.fits', outcube[*,*,nthis]
       STOP
    endif
    nthis = nthis + 1 
 
    if (debug gt 0) then STOP

  endfor                         ; end loop over frames in this cube (nim)

  endif else begin ;   if (nref_accept gt 0) and (nim ge njack) then begin
     print, '--------------' 
     print, '***WARNING***' 
     print, "Cube number: " + string(ic+1) + " no reference stars detected or "
     print, 'Too few frames in this cube:. There are ' + strn(nim) + ' frames, but'
     print, strn(njack) + ' are required for jackknife sampling.'
     print, 'Skipping cube ' + list[ic] + '.'
     print, "--------------" 
;     STOP
  endelse
  
  ; compute average and sigma
  ; require at least 3 valid frames per cube, otherwise
  ; the estimate of sigma will be very unreliable
  if (nthis gt 2) then begin
  
   ; trim cube 
   outcube = outcube[*,*,0:nthis-1]
   outmasks = outmasks[*,*,0:nthis-1]
   for x = 0, nax1-1 do begin
     for y = 0, nax2-1 do begin
      vals = reform(outcube[x,y,*])
     ; require at least 3 measurements per pixel
     ; for a reasonable estimation of sigma
      good = where(reform(outmasks[x,y,*]) gt 0, complement = bad, ngood)
      if ngood gt 2 then begin 
       vals = vals[good]
       wt[x,y] = wt[x,y] + n_elements(good)
      endif else begin
       vals[*] = 0
       wt[x,y] = 0
      endelse
      if ROBUST then begin
        RESISTANT_Mean, vals, Sigma_CUT, vals_mean, vals_sigma, Num_RejECTED
      endif else begin
        vals_mean = avg(vals)
        vals_sigma = stddev(vals)/sqrt(n_elements(vals)-1)
      endelse
      ssa[x,y] = vals_mean
      ssa_sigma[x,y] = vals_sigma
     endfor
    endfor
    writefits, tmpdir + 'ssa_' + strn(ic+1) + '.fits', ssa
    writefits, tmpdir + 'wt_' + strn(ic+1) + '.fits', wt
    writefits, tmpdir + 'ssa_sigma' + strn(ic+1) + '.fits', ssa_sigma
    name = list[ic]
    print, 'Averaged cube: ' + name
    wt_all = wt_all + wt
    ssa_all = ssa_all + wt*ssa
    ssa_sigma_all = ssa_sigma_all + wt*ssa_sigma^2
 
   ; Create jackknife samples
   ; -------------------------
   for j_c = 0, nthis-1 do begin
     i_j = (j_c+nt) mod njack
     ssa_jack[*,*,i_j] = ssa_jack[*,*,i_j] + outcube[*,*,j_c]
     wt_jack[*,*,i_j] = wt_jack[*,*,i_j] + outmasks[*,*,j_c]
   endfor
   nt = nt + nthis  

 endif  ; if (nthis gt 0)
 ic = ic + 1

endfor ; end loop over input cubes

good = where(wt_all gt 0, complement=bad)
ssa_all[good] = ssa_all[good]/wt_all[good]
if (bad[0] gt -1) then ssa_all[bad] = 0
ssa_sigma_all[good] = ssa_sigma_all[good]/wt_all[good]
if (bad[0] gt -1) then ssa_sigma_all[bad] = 0
writefits, outnam + '_wt.fits', wt_all, /COMPRESS
writefits, outnam + '.fits', ssa_all, /COMPRESS
writefits, outnam + '_sigma.fits', sqrt(ssa_sigma_all), /COMPRESS

; create jackknife images
; -----------------------
for i_j = 1, njack do begin
  wt = wt_jack[*,*,i_j-1]
  ssa = ssa_jack[*,*,i_j-1]
  good = where(wt gt 0)
  ssa[good] = ssa[good]/wt[good]
  writefits, outnam + '_s' + strn(i_j) + '.fits', ssa, /COMPRESS
endfor

END
 
