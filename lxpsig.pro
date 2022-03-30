PRO lxpsig, cube, maskcube, lxp, lxp_sigma, n_frames

Sigma_CUT = 3.0

sz = size(cube)
n1 = sz[1]
n2 = sz[2]
n_frames = sz[3]
lxp = fltarr(n1,n2)
lxp_sigma = fltarr(n1,n2)
for x = 0, n1-1 do begin
 for y = 0, n2-1 do begin
  mask = maskcube[x,y,*]
  valid = where(mask gt 0)
  if (valid[0] gt -1) then begin 
   vals = reform(cube[x,y,valid])
   RESISTANT_Mean, vals, Sigma_CUT, vals_mean, vals_sigma, Num_RejECTED
   lxp[x,y] = vals_mean
   lxp_sigma[x,y] = vals_sigma
  endif else begin
   lxp[x,y] = 0
   lxp_sigma[x,y] = 0
  endelse
 endfor
endfor

END
