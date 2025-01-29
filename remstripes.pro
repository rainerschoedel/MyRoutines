PRO REMSTRIPES, im, cleanim, hw, low_rej, high_rej, threshold, DEBUG = debug 

if not keyword_set(hw) then hw = 40
if hw lt 0 then hw = 40
if not keyword_set(low_rej) then low_rej = 10
if low_rej lt 0  then low_rej = 10
if not keyword_set(high_rej) then high_rej = 420
if high_rej lt 0 then high_rej = 420
if not keyword_set(threshold) then threshold = 5.0
if threshold lt 0  then threshold = 5.0

imcol_med, im, 0, low_rej, high_rej, outdata;, /trim    

; smooth median collapsed rows
smdata = medsmooth(outdata, hw)

; determine difference between median and smoothed median

ndata = size(outdata, /DIMENSIONS)
n = ndata[0]
corrdata = outdata - smdata

; Difference between median and smoothed median array
; allows to determine if counts are above threshold
; (ie are NOT 50Hz noise)

for i = 0, n - 1 do begin
    if abs(outdata[i] - smdata[i]) gt threshold then begin  ;
        corrdata[i] = 0.                                    ;
    endif
endfor
                   

; subtract 1d data column with 50Hz noise from every column of image

image_1d_sub, im, corrdata, 0, cleanim

if debug then begin
 sz = size(im)
 hhh1 = fltarr(sz[1],sz[2])
 hhh2 = fltarr(sz[1],sz[2])
 for x = 0, sz[1] -1 do begin 
  hhh1[x,*] = outdata
 endfor
 writefits, 'collapsed.fits', hhh1
 for x = 0, sz[1] -1 do begin 
  hhh2[x,*] = smdata
 endfor
 writefits, 'smooth.fits', hhh2
 writefits, 'diff.fits', hhh1 - hhh2
 writefits, 'im.fits', im
 writefits, 'cleanim.fits', cleanim
 STOP
endif

END
