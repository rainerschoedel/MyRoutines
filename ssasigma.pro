PRO SSASIGMA, indir, MASKDIR=maskdir, MASKLIST=masklist, outdir, OUTMASKDIR=outmaskdir, inlist, outim,  maskrad,  refsource, rebfac, debug, select, FWHM_SMOOTH=fwhm_smooth, OUTPUT = output, COMPRESS=compress, TMPDIR = tmpdir, ROBUST = robust, CEN_SHIFT = cen_shift, correct_sky = correct_sky, n_jack = n_jack

; PURPOSE: Make simple shift-and-add images from input image cubes
;          Return weight maps nad sigma maps.
; 
; LIMITAIONS: At least 3 valid frames per cube and 3 measurements per
;             pixel in the final images are required.
;             Only an integer pixel shift will be performed, even 
;             in the case of centroid shift.

; ROBUT:      If > 0 then use robust statistics for computing 
;             images of mean and standard deviation
;            THis approach will frequently produce NaN values! 
;
; CEN_SHIFT:  If > 0, shift to centroid, not to maximum pixel
;
; CORRECT_SKY: If >1 then a single constant is determined to be subtracted as sky offset from each frame

; VERSION DATE: 16 Dev 2015
;
; LAST CHANGES: 17 Jan 2025, Rainer Schoedel, added correct_sky keyword
;               21 Jan 2025, Rainer Schoedel, added n_jack jeyword

; Parameters for RESISTANT_Mean
Sigma_CUT = 3.0



withmask = keyword_set(maskdir)
if (not keyword_set(minsupp)) then minsupp = 0
if (not keyword_set(robust)) then robust = 0
if (not keyword_set(cen_shift)) then cen_shift = 0
if (not keyword_set(n_jack)) then n_jack = 0

name = ''
mname = ''
openr, inp, inlist, /get_lun
if withmask then openr, minp, masklist, /get_lun
ic = 0L

; determine size of image
; ------------------------
readf, inp, name
point_lun, inp, 0
cube = readfits(indir + name, NSLICE=1)
sz = size(cube)
nax1 = sz[1]
nax2 = sz[2]

; multiply with rebin factor
; --------------------------
maskrad = rebfac*maskrad
refsource = rebfac*refsource
nax1 = rebfac*nax1
nax2 = rebfac*nax2

; Define size and center of PSF
; -----------------------------
psfbox = (2*maskrad +20) + 1
cen = psfbox/2


; initialize variables
; ---------------------
ssa_all = fltarr(nax1,nax2)
ssa_sigma_all = fltarr(nax1,nax2)
wt_all = fltarr(nax1,nax2)
jack = fltarr(nax1,nax2,n_jack)
jack_wt = fltarr(nax1,nax2,n_jack)
jack_sigma = fltarr(nax1,nax2,n_jack)
nt = 0L
if (select eq 0) then $
  openw, out, 'stat.txt', /get_lun

while (not EOF(inp)) do begin

  readf, inp, name
  if withmask then readf, minp, mname

  ; read speckle data cube
  ; --------------------------
  cube = readfits(indir + name, header)
  sz = size(cube)
  if (sz[0] gt 2) then  nim = sz[3] else nim = 1
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  print, 'Read in cube: ' + name
  print, 'Number of frames in this cube: ' + strtrim(string(nim),2)
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  if withmask then maskcube = readfits(maskdir +  mname)
  sz = size(maskcube)
  if (sz[0] gt 2) then  nmask = sz[3] else nmask = 1
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  print, 'Read in cube: ' + name
  print, 'Number of frames in this cube: ' + strtrim(string(nim),2)
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  if (withmask and nmask eq 1) then begin
   mask = maskcube
   maskcube = fltarr(sz[1],sz[2],nim)
   for j = 0, nim-1 do maskcube[*,*,j] = mask
  endif
  ; if no mask is provided, create a full mask for the entire field
  ; it is later needed for weighting purposes
  if (not withmask) then begin
   maskcube = replicate(1,nax1,nax2,nim)
  endif

  nthis = 0L

  ; initialize images
  ssa = fltarr(nax1,nax2)
  wt = fltarr(nax1,nax2)
  ssa_sigma = fltarr(nax1,nax2)

  outcube = fltarr(nax1,nax2,nim)  
  outmasks = fltarr(nax1,nax2,nim)  

  for i = 0L, nim-1 do begin  ; start loop over frames in this cube

    ; load ith image 
    im = cube[*,*,i]
    mask = maskcube[*,*,i] 

    ; subtract constant sky?
    if correct_sky then begin
      sky_ind = where(mask gt 0)
      mmm, im(sky_ind), sky_mod, sky_sigma , sky_skew
      im = (im - sky_mod)*mask
    endif

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

; PSF extraction
; ################

    x_psf = refsource[0]
    y_psf = refsource[1]
    nref = n_elements(x_psf)
    xint = long(x_psf)
    yint = long(y_psf)
    psf = sub_array(im,psfbox,psfbox,REFERENCE=[xint,yint])
    psf = circ_mask(psf,cen,cen,maskrad,BORDER=0)
    if KEYWORD_SET(fwhm_smooth) then psf = filter_image(psf,FWHM_GAUSSIAN=fwhm_smooth)
    psfmax = max(psf, max_subscript)
    print, string(i+1) + ' ' + string(psfmax)
    if (select eq 0) then printf, out, psfmax
    mind = array_indices(psf,max_subscript)
    if (debug gt 0) then begin
      writefits, tmpdir + 'psf_preshift.fits', psf
    endif    

; shift each image and mask to position defined by PSF maximum or centroid
; ########################################################################

    if (psfmax gt select) then begin
;      if (debug gt 0) then begin
;       writefits, tmpdir + 'im.fits', im
;      endif
     ;centroid determination may go nuts if there are strong negativities
     ; set negative values = 0 to avoid this bug
     ; of course, this assumes that values< 0 are not important...
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
   endif ; if (psfmax gt select) 
  endfor  ; end loop over frames in this cube
  
 
  ; compute average and sigma
  ; require at least 3 valid frames per cube, otherwise
  ; the estimate of sigma will be very unreliable
  if (nthis gt 2) then begin
  
   ; trim cube (necessary if select > 0)
   if debug then begin
     writefits, tmpdir, 'outcube.fits', outcube
     writefits, tmpdir, 'outmasks.fits', outmasks
   endif
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
       if ROBUST then begin
         RESISTANT_Mean, vals, Sigma_CUT, vals_mean, vals_sigma, Num_RejECTED
       endif else begin
         vals_mean = avg(vals)
         vals_sigma = stddev(vals)/sqrt(ngood)
       endelse
         ssa[x,y] = vals_mean
         ssa_sigma[x,y] = vals_sigma
         wt[x,y] = wt[x,y] + ngood
       endif else begin
         ssa[x,y] = 0
         ssa_sigma[x,y] = 0
       endelse
      endfor
    endfor
    writefits, tmpdir + 'ssa_' + strn(ic+1) + '.fits', ssa
    writefits, tmpdir + 'wt_' + strn(ic+1) + '.fits', wt
    writefits, tmpdir + 'ssa_sigma' + strn(ic+1) + '.fits', ssa_sigma
    print, 'Averaged cube: ' + name
    wt_all = wt_all + wt
    ssa_all = ssa_all + wt*ssa
    ssa_sigma_all = ssa_sigma_all + wt*ssa_sigma^2
 
    ; create jackknife samples
    skip = (ic mod n_jack) ; "ic" is subscript of  current cube
    for i_j = 0, n_jack-1 do begin
      if (i_j ne skip) then begin ; and (i lt n_jack) then begin
            jack[*,*,i_j] = jack[*,*,i_j] +  wt*ssa
            jack_wt[*,*,i_j] = jack_wt[*,*,i_j] + wt
            jack_sigma[*,*,i_j] = jack_sigma[*,*,i_j] + wt*ssa_sigma^2
            print, 'Cube ' + strn(ic) + ' to jack sample ' + strn(i_j)
      endif else print, 'Skipped ' + strn(ic)
    endfor
   if (output gt 0) then begin 
    writefits, outdir + name, outcube[*,*,0:nthis-1], header
    writefits, outmaskdir + mname , outmasks[*,*,0:nthis-1], COMPRESS=compress
   endif
;writefits, tmpdir + 'wt_jack.fits', jack_wt
;writefits, tmpdir + 'jack.fits', jack

 endif  ; if (nthis gt 0)
 ic = ic + 1

endwhile ; end loop over input cubes

if withmask then free_lun, minp
if (select eq 0) then free_lun, out
free_lun, inp

good = where(wt_all gt 0, complement=bad)
ssa_all[good] = ssa_all[good]/wt_all[good]
if (bad[0] gt -1) then ssa_all[bad] = 0
ssa_sigma_all[good] = ssa_sigma_all[good]/wt_all[good]
if (bad[0] gt -1) then ssa_sigma_all[bad] = 0
writefits, outim + '_wt.fits', wt_all, /COMPRESS
writefits, outim + '.fits', ssa_all, /COMPRESS
writefits, outim + '_sigma.fits', sqrt(ssa_sigma_all), /COMPRESS

; create jackknife images

for i_j = 0, n_jack-1 do begin
  jack_im_wt = jack_wt[*,*,i_j]
  good = where(jack_im_wt gt 0)
  jack_im = jack[*,*,i_j]
  jack_im[good] = jack_im[good]/jack_im_wt[good]
  jack_noise = jack_sigma[*,*,i_j]
  jack_noise[good] = jack_noise[good]/jack_im_wt[good]
  jack_noise = sqrt(jack_noise)
  writefits, outim + '_s' + strn(i_j+1) + '.fits', jack_im, /COMPRESS
  writefits, outim + '_s' + strn(i_j+1) + ' _sigma.fits', jack_noise, /COMPRESS
 endfor

END

