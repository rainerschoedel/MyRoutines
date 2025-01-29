PRO REMOVE_NEARBY, x, y, f, mag, x_psf, y_psf, mag_psf, delta_mag, mrad

; PURPOSE
; De-select possible spurious sources in diffraction spikes of bright
; stars from the list of detected stars so that they will be used in
; PSF estimation
  
  
; INPUT
; X, Y, MAG: XY positions and magnitudes of all stars
; X_ISO, Y_ISO, MAG_ISO"XY positions and magnitudes of
; PSF stars
; 
; delta_mag: magnitude difference criterion
; delta_r: distance criterion 

; Only for test purposes. Adapt to image size.
;  n1 = 2048
;  n2 = 2048

 print, 'Remove nearby stars.'
 print, 'Input ' + strn(n_elements(x_psf)) + ' PSF stars. '
 print, 'Input ' + strn(n_elements(x)) + ' stars. '
 dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 1.0, Sigma_y: 1.0, Angle: 0.0})
 n_psf = n_elements(mag_psf)
 ; now check deselection criteria
 for s = 0, n_psf-1 do begin
  ind_keep = indgen(n_elements(x))
  d = sqrt((x_psf[s] - x)^2 + (y_psf[s] - y)^2)
  nearby = where(d gt mrad/10 and d le mrad, n_near)
  if (n_near gt 0) then begin
     diff_mag = mag[nearby] - mag_psf[s]
     deselect = where(diff_mag gt delta_mag, count)
     if (count gt 0) then begin
        REMOVE, nearby[deselect], ind_keep
        x = x[ind_keep]
        y = y[ind_keep]
        f = f[ind_keep]
        mag = mag[ind_keep]
;        im = image_model(x,y,f,n1,n2,'gaussian', dat)
;        writefits, 'tmp/map_' + strn(s) + '_cleaned.fits', im, /COMPRESS
     endif
   endif
endfor
 
 n_keep = n_elements(ind_keep)
 print, 'Keep ' + strn(n_keep) + ' stars. '

END
