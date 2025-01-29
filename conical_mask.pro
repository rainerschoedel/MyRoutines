FUNCTION CONICAL_MASK, array, x0, y0, tau0, w_tau

; NAME:
;	CONICAL_MASK
;
; PURPOSE:
;	Apply a conical mask to a 2D array, setting to a pre-fixed value all
;	the pixels that lie within a pre-defined cone of a reference
;	pixel (apex) and to zero all pixels outside that cone.
;
; CATEGORY:
;	Signal processing.
;
; CALLING SEQUENCE:
;	Result = CONICAL_MASK(Array, X0, Y0, Theta0, W_THETA)
;
; INPUTS:
;	Array:	2D array to mask
;
;	X0, Y0:	Coordinates of apex of cone
;
;	Theta0:	Direction of cone in degrees, increasing
;	counter-clockwise and zero toward positive y-axis
;  
;       W_Theta: Width of cone
;
;
; OUTPUTS:
;	Result:	Array with region defined by conical mask set to 1 
;               and region outside that mask set to zero
;
; MODIFICATION HISTORY:
; 	Written by:	Rainer Schoedel, August 2013.
;
;       edited if statement and changed to 
;       "if (keyword_set(bord) eq0)"
;       Rainer Schoedel October 2010


	on_error, 2
	sz = size(array)
        n1 = sz[1]
        n2 = sz[2]


        tauhi = tau0 + w_tau/2.
        taulo = tau0 - w_tau/2.
        ; special cases: talo <0 or tauhi > 360.
        if (taulo lt 0) then taulo_mod = 360. + taulo
        if (tauhi gt 360.) then tauhi_mod = tauhi - 360.
 
        cone = array
        cone[*,*] = 0
        for ix = 0L, n1-1 do begin
         for iy = 0L, n2-1 do begin
          dy = iy - y0
          dx = ix - x0
          tau = 180./!PI * atan(dy,dx)
          ; convert to angle definition 0 <= tau <= 360 
          ; tau=0 toward positive y
          ; increasing counter-clockwise
          if (tau lt 0) then tau = tau + 360. ; define from 0 o 360 degrees
          tau = tau - 90.                     ; positive toward positive y
          if (tau lt 0) then tau = tau + 360. ; only positive angles
          if (tau le tauhi and tau ge taulo) then cone[ix,iy] = 1
          ; special cases: talo <0 or tauhi > 360.
          if (taulo lt 0) then begin
           if (tau le tauhi or tau ge taulo_mod) then cone[ix,iy] = 1
          endif
          if (tauhi gt 360.) then begin
           if (tau le tauhi_mod or tau ge taulo) then cone[ix,iy] = 1
          endif
         endfor
        endfor

	return, cone * array

end



END
