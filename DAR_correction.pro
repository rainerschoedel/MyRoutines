pro DAR_apply, meteo, frames

; Original code writen by Kora Muzic and modified by Rainer Schoedel 
; in February 2011.
; 
; INPUT:
; meteo:  file from ESO data base http://archive.eso.org/wdb/wdb/eso/ambient_paranal/form
;         with meteorological information
; frames: file with list of FITS - images for which the DAR correction wil be
;         calculated; FITS headers must contain UT and LST


set_plot,'ps',/interpolate
device, XSIZE=15, YSIZE=15, XOFFSET=0, YOFFSET=0, $
        FILENAME='DAR2D.eps',$
        /portrait, /color,BITS_PER_PIXEL=8, encapsulated=0

!P.MULTI=[0,0,0]
!P.CHARSIZE=1.3
!P.THICK=2.0
!P.CHARTHICK=2.0

;use this procedure to apply DAR to the Hubble data, in
;order to simulate how they would look like if we were 
;observing from the ground

;Paranal UT4 latitude
obslat=ten(-24,37,31)

;observing wavelength in microns
lambda=2.18

readcol, meteo,night,int,T30,T2,Tg,DT30,DT2,RH30,RH2,P2,delimiter='|',format='(A,F,F,F,F,A,F,A,F,F)'
readcol, frames,framelist,format='A'

for i=0,n_elements(framelist)-1 do begin
;for i=0, 0 do begin  
  header=HEADFITS(framelist[i])

 ;extract UT and elevation angle from the frame header
  time=SXPAR(header,'UT')
  LST=SXPAR(header,'LST') ;seconds!
; elevation=double(esosxpar(header,'TEL ALT'))
  
 ;reference position: guide star
  RAref=ten(00,24,11.0)*15.0
  DECref=ten(-72,05,52.0)

 ;convert LST to degrees
  LST=LST/3600.0*15.0

 ;read the atmospheric conditions at the given UT
  aa = float(STRSPLIT(time, ':',/REGEX, /EXTRACT))
  time= ten(aa[0],aa[1],aa[2])
 
  tt=dblarr(n_elements(night))
  for j=0, n_elements(night)-1 do begin
   aa = STRSPLIT(night[j], ' ',/REGEX, /EXTRACT)
   aa = STRSPLIT(aa[1], ':',/REGEX, /EXTRACT)
   tt[j]=ten(aa[0],aa[1],aa[2])
  endfor

  tmp=replicate(time, n_elements(night))
  tmp=abs(tmp-tt)
  tmpmin=min(tmp,sub)

  T=T2[sub]
  RH=RH2[sub]/100.0
  P=P2[sub]
 
 ;Hubble plate scale in arcsec/pixel 
  platescale=0.043 

 ;orientation of the Hubble detector (E of N)
  PosAng=0.0*!DTOR
  
  readcol, 'HST.positions',x,y
  xref=512.0
  yref=512.0
  
 ;first de-rotate (N=y, E=-x), around the reference point
  xderot=xref+(x-xref)*cos(-PosAng)-(y-yref)*sin(-PosAng)
  yderot=yref+(x-xref)*sin(-PosAng)+(y-yref)*cos(-PosAng)

; convert x, y coordinates to celestial coordinates
; needs to be re-written
  RA = RAref - (xderot-xref)*platescale/3600.0/cos(DECref*!DTOR); ??????
  Dec= DECref+ (yderot-yref)*platescale/3600.0

 ;hour angle 
  HA = LST - RA ;in degrees

 ;parallactic angle (angle between North direction and Zenith) 
  ParAng = atan(sin(HA*!DTOR)/(cos(dec*!DTOR) * tan(obslat*!DTOR) - sin(dec*!DTOR) * cos(HA*!DTOR)))

 ; zd=acos(sin(obslat*!DTOR)*sin(dec*!DTOR)+cos(obslat*!DTOR)*cos(dec*!DTOR)*cos(HA*!DTOR)))

 ;atmospheric refraction angle at the reference position
 ;DARref=DAR(obslat, LST, RAref, decref, P, T, RH, lambda)
 
 ; compute change of position of each star due to atmospheric 
 ; refraction
 ; DAR is given in arcseconds
 ; DAR=abs(DAR(obslat, LST, RA, dec, P, T, RH, lambda)-DARref)
  DAR=DAR(obslat, LST, RA, dec, P, T, RH, lambda)

  ; calculate offset on detctor due to AR
  xDAR=DAR*sin(PosAng-ParAng)/platescale
  yDAR=DAR*cos(PosAng-ParAng)/platescale
  
  ; appply the atmospheric refraction to the HST pixel coordinates
  xnew=x+xDAR
  ynew=y+yDAR

  aa = STRSPLIT(framelist[i], '/',/REGEX, /EXTRACT)
  openw, l1, 'HST_DAR_'+aa[n_elements(aa)-1],/get_lun
  for j=0,n_elements(x)-1 do printf, l1,xnew[j],ynew[j]
  free_lun, l1

  plot_2D, PosAng,xref,yref,ParAng,xnew,ynew,DAR
endfor

device,/close

end

pro plot_2D, PosAng, xref, yref,ParAng,x,y,DAR

 
  readcol,'HST.positions',x0,y0
 ; readcol,'HST_DAR_im27.fits',x,y

  plot,findgen(1024),findgen(1024),/nodata,/isotropic,xrange=[0,1024],yrange=[0,1024]

  plots, xref, yref,psym=1,color=fsc_color('red')

  sc=10
  for j=0,n_elements(x0)-1 do begin
    ARROW, x0[j],y0[j],(x[j]+sc*x0[j])/(1+sc),(y[j]+sc*y0[j])/(1+sc),hsize=180.0 ,/data
    xyouts, x0[j],y0[j],strtrim(DAR[j],2)
  endfor

  x=findgen(2000)-1000
  k=tan(!PI/2+PosAng+ParAng[0])
  oplot,x,k*x+yref-k*xref,linestyle=1,color=fsc_color('dark grey')

end

function DAR, lat, LST, RA, dec, p, T, H, lambda

;standard values
Ts=288.15
ps=1013.25 ;hPa

;Maximum water vapor pressure as a function of the temperature (table
;1 from Helminiak 2009
;Here Temperature has to be in degrees celsius!
;Formula/interpolation see Section 3 in Helminiak 2009 
 t0=[-45,   -40,   -35, -30,-25,-20, -15,-10, -5,  0,    5,   10,  15,   20,   25,    30,  35,  40,   45,    50]
pwmax=[0.108,0.185,0.309,0.5,0.8,1.25,1.9,2.68,4.21,6.11,8.72,12.28,17.05,23.27,31.66,42.41,56.2,73.72,95.77,123.3]
pw=H*spline(t0,pwmax,T)


;refractive index of air
n=1+(64.328+29498.1/(146-lambda^(-2.0D))+255.4/(41-lambda^(-2.0D)))*(p*Ts/ps/(T+273.15))*10^(-6.0D) $
     - 43.49*(1-0.007956/lambda^2.0D)*(pw/ps)*10^(-6.0)

;true zenith distance
;obslat=observer's latitude
;targdec=declination of the target
;HA=hour angle of the target
;LST=local sidereal time

HA = LST - RA ;in degrees

; calculate true zenith distance: (4) in Roe 2002
zt=acos(sin(lat*!DTOR)*sin(dec*!DTOR)+cos(lat*!DTOR)*cos(dec*!DTOR)*cos(HA*!DTOR))

;result in arcsec: (3) in Roe 2002
DAR=206265.0D*((n^2-1)/(2.0*n^2))*tan(zt)

return, DAR
end





;"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

pro DAR_plot

;this is to produce the plot that shows the amplitude of DAR vs
;elevation, for several source separations
;It reproduces the Figure 4 from Yelda et al. (2010)
;calls the procedure DAR_elev.pro which is similar to DAR, but simpler
;because you just supply an array of elevations (no RA, DEC, LST etc.)

set_plot,'ps',/interpolate
device, XSIZE=15, YSIZE=15, XOFFSET=0, YOFFSET=0, $
        FILENAME='DAR.eps',$
        /portrait, /color,BITS_PER_PIXEL=8, encapsulated=0
!P.MULTI=[0,0,0]
!P.CHARSIZE=1.3
!P.THICK=2.0
!P.CHARTHICK=2.0
plot, findgen(90),findgen(90),/nodata,yrange=[0,8],xrange=[20,90],xtitle='elevation [deg]',ytitle='DAR [mas]',$
      XTHICK=2.0, YTHICK=2.0,xstyle=1,ystyle=1

;Paranal UT4 latitude
lat=ten(-24,37,31)

;observing wavelength in microns
lambda=2.18

elev=findgen(90)+1
 
T=9.38
H=0.1105
P=742.33

dist=[1,5,10] ;arcsec

for i=0,n_elements(dist)-1 do begin
 DAR1=DAR_elev(elev,p,T,H,lambda)
 
 elev2=double(elev+dist[i]/3600.0)
 DAR2=DAR_elev(elev2,p,T,H,lambda)
 
 DR=abs(DAR1-DAR2)*1000.0
 oplot,elev,DR
 xyouts, 70, DR(where(elev eq 70))+0.03, 'd='+strtrim(string(dist(i)),2)+'"'
 print, dist[i],DR(where(elev eq 40))
endfor


;check how the absolute value of the
;refraction angle depends on zenith distance
 DAR=DAR_elev(elev,p,T,H,lambda)
 
 plot, 90-elev, DAR, xtitle='zenith distance [deg]',ytitle='atmospheric refraction [acsec]'

device,/close
end



function DAR_elev, elev, p, T, H, lambda

;standard values
Ts=288.15
ps=1013.25 ;hPa

;Maximum water vapor pressure as a function of the temperature (table
;1 from Helminiak (2009), New Astronomy 14 521-527
t0=[-45,   -40,   -35, -30,-25,-20, -15,-10, -5,  0,    5,   10,  15,   20,   25,    30,  35,  40,   45,    50]
pwmax=[0.108,0.185,0.309,0.5,0.8,1.25,1.9,2.68,4.21,6.11,8.72,12.28,17.05,23.27,31.66,42.41,56.2,73.72,95.77,123.3]
pw=H*spline(t0,pwmax,T)

;refractive index of air
n=1+(64.328+29498.1/(146-lambda^(-2.0D))+255.4/(41-lambda^(-2.0D)))*(p*Ts/ps/(T+273.15))*10^(-6.0D) - 43.49*(1-0.007956/lambda^2.0D)*(pw/ps)*10^(-6.0)

zt=90.0-elev
;result in arcsec
DAR=206265.D*((n^2-1)/2.0/n^2)*tan(zt*!DTOR)
return, DAR

end



pro DAR_correct
;this should be a procedure to correct for the DAR effects
;It's not finished!!!

;Paranal UT4 latitude
obslat=ten(-24,37,31)

;observing wavelength in microns
lambda=2.18

readcol, 'meteo_20100825.dat',night,int,T30,T2,Tg,DT30,DT2,RH30,RH2,P2,delimiter='|',format='(A,F,F,F,F,F,F,F,F,F)'
readcol, 'framelist.dat',framelist,format='A'

for i=0,n_elements(framelist)-1 do begin
;for i=0, 0 do begin  
  header=HEADFITS(framelist[i])

 ;extract UT and elevation angle from the frame header
  time=SXPAR(header,'UT')
  LST=SXPAR(header,'LST') ;seconds!
  
 ;reference position: guide star
  RAref=ten(00,24,20.16)*15.0
  DECref=ten(-72,05,38.0)

 ;convert LST to degrees
  LST=LST/3600.0*15.0

 ;read the atmospheric conditions at the given UT
  aa = float(STRSPLIT(time, ':',/REGEX, /EXTRACT))
  time= ten(aa[0],aa[1],aa[2])
 
  tt=dblarr(n_elements(night))
  for j=0, n_elements(night)-1 do begin
   aa = STRSPLIT(night[j], ' ',/REGEX, /EXTRACT)
   aa = STRSPLIT(aa[1], ':',/REGEX, /EXTRACT)
   tt[j]=ten(aa[0],aa[1],aa[2])
  endfor

  tmp=replicate(time, n_elements(night))
  tmp=abs(tmp-tt)
  tmpmin=min(tmp,sub)

  T=T2[sub]
  RH=RH2[sub]/100.0
  P=P2[sub]
 
 ;Hubble plate scale in arcsec/pixel 
  platescale=0.043 
 
 ;hour angle 
  HA = LST - RA ;in degrees

 ;parallactic angle (angle between North direction and Zenith) 
  ParAng = atan(sin(HA*!DTOR)/(cos(dec*!DTOR) * tan(obslat*!DTOR) - sin(dec*!DTOR) * cos(HA*!DTOR)))

;calculate refractive index of air
 ;standard values
 Ts=288.15
 ps=1013.25 ;hPa

 ;Maximum water vapor pressure as a function of the temperature (table
 ;1 from Helminiak 2009
 ;Here Temperature has to be in degrees celsius!
 ;Formula/interpolation see Section 3 in Helminiak 2009 
 t0=[-45,   -40,   -35, -30,-25,-20, -15,-10, -5,  0,    5,   10,  15,   20,   25,    30,  35,  40,   45,    50]
 pwmax=[0.108,0.185,0.309,0.5,0.8,1.25,1.9,2.68,4.21,6.11,8.72,12.28,17.05,23.27,31.66,42.41,56.2,73.72,95.77,123.3]
 pw=H*spline(t0,pwmax,T)

 ;refractive index of air
 n=1+(64.328+29498.1/(146-lambda^(-2.0D))+255.4/(41-lambda^(-2.0D)))*(p*Ts/ps/(T+273.15))*10^(-6.0D) - 43.49*(1-0.007956/lambda^2.0D)*(pw/ps)*10^(-6.0)
 
 zt=acos(sin(obslat*!DTOR)*sin(dec*!DTOR)+cos(obslat*!DTOR)*cos(dec*!DTOR)*cos(HA*!DTOR))
 el=asin(sin(dec*!DTOR)*sin(obslat*!DTOR)+cos(obslat*!DTOR)*cos(dec*!DTOR)*cos(HA*!DTOR))





endfor
end


