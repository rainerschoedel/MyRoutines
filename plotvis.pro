PRO PLOTVIS, filename

; Read an OIFITS file and plot
; V2 over baseline
; V2 in UV plane

; read OIFITS file
read_oidata, filename, oiarray, oitarget, oiwavelength, oivis, oivis2, oit3

; get baselines
; u and v are in meters
u = oivis2.ucoord
v = oivis2.vcoord
r = sqrt(u^2 + v^2)
n_baselines = n_elements(r)

; get visibilities
v2 = fltarr(n_baselines)
v2err = fltarr(n_baselines)
for i = 0, n_baselines-1 do begin
  v2[i] = *oivis2[i].vis2data
  v2err[i] = *oivis2[i].vis2err
endfor

; other information
tname = oitarget[0].target

; V2 over baseline
ps_settings, 'v2r.eps', /eps, pmulti=[0,1,1], axisratio=3./4.,charsize=0.3
plot, r, v2, xran=[0.,8.], yran=[0.,1.],title=tname,xtitle='baseline length[m]',ytitle='V^2 calibrated', /nodata
oploterr, r, v2, v2err,4
ps_close

; V2 over baseline, coded for baseline angle
theta = atan(v,u)/!PI * 180.0
ix = where(theta lt 0)
theta[ix] = theta[ix] + 360.0
ps_settings, 'v2rtheta.eps', /eps, pmulti=[0,1,1], axisratio=3./4.,charsize=0.3
plot, theta, v2, xran=[0.,360.], yran=[0.,1.],title=tname,xtitle='baseline angle [deg]',ytitle='V^2 calibrated', /nodata
ix = sort(r)
sr = r[ix]
stheta = theta[ix]
sv2 = v2[ix]
sv2err = v2err[ix]
sr = 0.01*round(100*sr)
bl = uniq(sr)
nb = n_elements(bl)
print, bl
for j = 0L, nb-1 do begin
  if (j eq 0) then six = indgen(bl[0]+1) $
  else six = 1 + bl[j-1] + indgen(bl[j]-bl[j-1])
  sym = j+1
  loadct, 34
  oplot, stheta[six], sv2[six],psym=4,color=(20*(j+1))
  xyouts, 20, 0.05*(1+j), string(sr[bl[j]]),color=(20*(j+1))
endfor
ps_close


END
