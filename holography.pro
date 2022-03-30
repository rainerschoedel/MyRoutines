PRO CALEVAL, X, P, YMOD

  YMOD = P[0] + P[1]*X
  YMOD = image_shift(YMOD,P[2],P[3])

END


;---------------------

PRO HOLOGRAPHY,  indir, list_numer, list_masknumer, list_denom, OUTNAM = outnam, AIRY=airy, OUT_ITER=out_iter, tmpdir=tmpdir, NSUB=nsub, MINSUPP=minsupp, MAXSUPP = maxsupp, RAWOUT = rawout, DEBUG = debug

; Created on August 11, 2017
; HOLO_MOSIAC must be run first to create the necessary input files

if (not KEYWORD_SET(debug)) then debug = 0
if (not KEYWORD_SET(maxsupp)) then maxsupp = 1.1
if (not KEYWORD_SET(minsupp)) then minsupp = 0.01
if (not KEYWORD_SET(rawout)) then rawout= 0
if (not KEYWORD_SET(smoothmask)) then smoothmask = 0

n_cubes = n_elements(list_numer)

restore, indir + list_numer[0]
sz = size(ghq_arr)
nax1 = sz[1]
nax2 = sz[2]

cen1 = nax1/2
cen2 = nax2/2
center = [cen1,cen2]

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
fimasknumer = fltarr(nax1,nax2)
fimasknumersubreal = fltarr(nax1,nax2,nsub)
fimasknumersubim = fltarr(nax1,nax2,nsub)

for ic = 0, n_cubes-1 do begin  ; start loop over all input cubes

  ; initialize variables for holography
  ;---------------------------------------

  mask = replicate(0.0,nax1,nax2)
;  mask_sub = replicate(0.0,nax1,nax2,nsub)
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

  ; read data cubes
  ; --------------------------

   restore, indir + list_numer[0]
   restore, indir + list_denom[0]
   restore, indir + list_masknumer[0]
   sz = size(ghq_arr)
   nim = sz[3]
   
   for i = 0L, nim-1 do begin  ; start loop over frames in this cube

      GHQ = ghq_arr[*,*,i]
      GMHQ = gmhq_arr[*,*,i]
      H2 = h2_arr[*,*,i]
      masknumer = masknumer + GMHQ
      fimasknumer = fimasknumer + GMHQ
    
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
      masknumersubreal[*,*,n123] = masknumersubreal[*,*,n123] + REAL_PART(GMHQ)
      masknumersubim[*,*,n123] = masknumersubim[*,*,n123] + IMAGINARY(GMHQ)
        
      n123 = ngood mod nsub
      finumersubreal[*,*,n123] = finumersubreal[*,*,n123] + REAL_PART(GHQ)
      fidenomsubreal[*,*,n123] = fidenomsubreal[*,*,n123] + REAL_PART(H2)
      finumersubim[*,*,n123] = finumersubim[*,*,n123] + IMAGINARY(GHQ)
      fidenomsubim[*,*,n123] = fidenomsubim[*,*,n123] + IMAGINARY(H2)
      fimasknumersubreal[*,*,n123] = fimasknumersubreal[*,*,n123] + REAL_PART(GMHQ)
      fimasknumersubim[*,*,n123] = fimasknumersubim[*,*,n123] + IMAGINARY(GMHQ)
  
      hhh = REAL_PART(FFT(OTF*(GHQ/H2), /INVERSE))
      rcim = shift(hhh, center+1)

      if (debug gt 0) then begin
         mmm = REAL_PART(FFT(OTF*(GMHQ/H2), /INVERSE))
         rcmask = shift(mmm, center+1)
         writefits, tmpdir + 'rcmask.fits', rcmask
         writefits, tmpdir + 'rcim.fits', rcim
       endif
  
      print, "Frame number: " + string(i+1) + " ok"

      if (debug gt 0) then STOP

    ; output of control data
    ; ----------------------

    if ((nthis mod out_iter) eq 0 and nthis gt 0) then begin
     hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
     rcim = shift(hhh, center+1)
     mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
     rcmask = shift(mmm, center+1)
     writefits, tmpdir + 'rcmask.fits', rcmask
     support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
     rcim[support] = rcim[support]/rcmask[support]
     if (nosup[0] gt -1) then begin
      rcim[nosup] = 0
     endif
     writefits, tmpdir + 'rcim_tmp.fits', rcim
    endif

  endfor                         ; end loop over frames in this cube

  ; output of holography image for this cube
  ; ----------------------------------------------------
  hhh = REAL_PART(FFT(OTF*(numer/denom), /INVERSE))
  rcim = shift(hhh, center+1)
  mmm = REAL_PART(FFT(OTF*(masknumer/denom), /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
  if (support[0] gt -1) then begin
   rcim[support] = rcim[support]/rcmask[support]
   if (nosup[0] gt -1) then begin
    rcim[nosup] = 0
    rcmask[nosup] = 0
   endif
   writefits, tmpdir + 'mask_' + strtrim(string(ic+1),2)+'.fits', mask
   writefits, tmpdir + 'support_' + strtrim(string(ic+1),2)+'.fits', rcmask
  endif
  writefits, tmpdir + 'holo_' + strtrim(string(ic+1),2)+'.fits', rcim, header

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
   writefits, tmpdir + 'support_' + strtrim(string(ic+1),2) + '_s' + strtrim(string(is),2) +'.fits', rcmask
   writefits, tmpdir + 'holo_' + strtrim(string(ic+1),2)+ '_s' + strtrim(string(is),2) + '.fits', rcim
  endfor

 ; output of current total holography image
  ; ----------------------------------------------------
   hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
  rcim = shift(hhh, center+1)
  mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
  rcim[support] = rcim[support]/rcmask[support]
  if (nosup[0] gt -1) then begin
   rcim[nosup] = 0
   rcmask[nosup] = 0
  endif
  writefits, tmpdir + 'current_mask.fits', fimask
  writefits, tmpdir + 'current_support.fits', rcmask
  writefits, tmpdir + 'current.fits', rcim

endfor ; end loop over input cubes


; final image
; --------------------
hhh = REAL_PART(FFT(OTF*(finumer/fidenom), /INVERSE))
rcim = shift(hhh, center+1)
mmm = REAL_PART(FFT(OTF*(fimasknumer/fidenom), /INVERSE))
rcmask = shift(mmm, center+1)
support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
rcim[support] = rcim[support]/rcmask[support]
if (nosup[0] gt -1) then begin
 rcim[nosup] = 0.0
 rcmask[nosup] = 0
endif
writefits, outnam+'_support.fits', rcmask
writefits, outnam + '_mask.fits', fimask
writefits, outnam+'.fits', rcim
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
 writefits, outnam+'_support_raw.fits', rcmask
 writefits, outnam + '_mask_raw.fits', fimask
 writefits, outnam+'_raw.fits', rcim
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
  maskhhh = complex(fimasknumersubreal[*,*,is-1],fimasknumersubim[*,*,is-1])
  mmm = REAL_PART(FFT(maskhhh/dhhh, /INVERSE))
  rcmask = shift(mmm, center+1)
  support = where(rcmask gt minsupp and rcmask lt maxsupp, complement=nosup)
  rcim[support] = rcim[support]/rcmask[support]
  if (nosup[0] gt -1) then begin
   rcim[nosup] = 0
   rcmask[nosup] = 0
  endif
  writefits, outnam + '_s_raw' + strtrim(string(is),2) + '.fits', rcim
 endfor
endif

END
 
