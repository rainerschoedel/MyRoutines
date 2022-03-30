PRO ISOLATED_STARS, x, y, mag, x_iso, y_iso, mag_iso, delta_mag, delta_r, ind_iso

; INPUT
; X, Y, MAG: XY positions and magnitudes of all stars
; X_ISO, Y_ISO, MAG_ISO"XY positions and magnitudes of
; candidate isolated stars
; 
; delta_mag: magnitude difference criterion
; delta_r: distance criterion 

 ind_iso = [] ; make sure that ind_iso is empty
 n_iso = n_elements(mag_iso)
 ; now check isolation criteria
 for s = 0, n_iso -1 do begin
  d = sqrt((x_iso[s] - x)^2 + (y_iso[s] - y)^2)
  nearby = where(d gt 0 and d le delta_r, n_near)
  if (n_near lt 1) then begin
   if n_elements(ind_iso) gt 0 then begin
    ind_iso = [ind_iso,s]
   endif else begin
    ind_iso = s
   endelse
  endif else begin
   diff_mag = mag[nearby] - mag_iso[s]
   if min(diff_mag) gt delta_mag then begin
    if n_elements(ind_iso) gt 0 then begin
     ind_iso = [ind_iso,s]
    endif else begin
     ind_iso = s
    endelse
   endif 
  endelse
 endfor
 
 n_iso = n_elements(ind_iso)
 print, 'Found ' + strn(n_iso) + ' isolated stars. '

END
