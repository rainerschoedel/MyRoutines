PRO CALEVAL, X, P, YMOD

  YMOD = P[0] + P[1]*X
  YMOD = image_shift(YMOD,P[2],P[3])

END


;---------------------

FUNCTION EXTRACT_PSF, im, x_psf, y_psf, boxhw, maskrad, normrad, NSIGMA=nsigma, REBFAC=rebfac, STARLIST=starlist, UNWEIGHTED = unweighted, DEBUG = debug, SUBPIX=subpix, ITER = iter

; VERSION DATE: 11 March 2014

; Define some defaut values
if  not keyword_set(nsigma) then nsigma = 3.0
if  not keyword_set(rebfac) then rebfac = 1.0
if  not keyword_set(unweighted) then unweighted = 0
if  not keyword_set(debug) then debug = 0
if  not keyword_set(subpix) then subpix = 1
if  not keyword_set(iter) then iter = 0

; number of iterations
 iter = n_elements(nsigma) - 1

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

; apply rebin factor
; use new variabels for rebinned quantities so that the parameters
; in the callin gprogram will not be changed
mrad = rebfac*maskrad
nax1 = rebfac*nax1
nax2 = rebfac*nax2
nrad = rebfac*normrad
bhw = rebfac*boxhw
psfim = FREBIN(im,nax1,nax2,/TOTAL)

; numbers for subarray extraction
boxsize = long(2*bhw + 1)
cen1 = nax1/2
cen2 = nax2/2
center = [cen1,cen2]

; dummy needed  ringmask for background estimation
dummy = fltarr(boxsize,boxsize)
dummy[*,*] = 1
bgring = circ_mask(dummy,bhw,bhw,mrad,BORDER=0, INNER=1)
bgind = where(bgring gt 0)

; reference sources
x_ref = rebfac * x_psf
y_ref = rebfac * y_psf
nref = n_elements(x_ref)
xint = round(x_ref)
yint = round(y_ref)

; PSF extraction
; ################

sub_arrays, psfim, xint, yint, boxsize, stack, masks
; save non-normalized stack for photometry
rawstack = stack

for iref = 0, nref-1 do begin
  subim = stack[*,*,iref]
  if (subpix gt 0) then begin
    subim = image_shift(subim,xint[iref]-x_ref[iref],yint[iref]-y_ref[iref])
  endif
  ; subtract local background and normalize sub-image 
  normim = subim - median(subim[bgind])
  normim = circ_mask(normim, bhw, bhw, nrad,BORD=0)
  if (unweighted eq 1 and total(normim) gt 0) then begin
    subim = subim/total(normim)
  endif
  stack[*,*,iref] = subim
 endfor
 psf = stack_avg(stack, MASK=masks, AVGYN=0)
 
 ; subtract remaining background, suppress noise,  and mask PSF
 ; -----------------------------------------------------------
   if (debug gt 0) then begin
     writefits, tmpdir + 'im.fits', psfim
     writefits, tmpdir + 'bgring_0.fits', psf*bgring
     writefits, tmpdir + 'psf_0raw.fits', psf
     writefits, tmpdir + 'stack_0.fits', stack*masks
    endif
    psf = psf - median(psf[bgind])  ; subtract background
    noise = stddev(psf[bgind])
    suppress = where(psf lt nsigma[0] * noise, count)        
    if (count gt 0) then psf[suppress] = 0
    psfcen = centroid(circ_mask(psf,bhw,bhw,mrad,BORDER=0))   
    psf = circ_mask(psf,psfcen[0],psfcen[1],mrad,BORDER=0)
    psf = psf/total(psf)
 
 return, psf

END
 
