PRO srebin,image,rmap,dmap,nsout,nlout
;+
; NAME:
;   SREBIN
;
; PURPOSE:
;   Shrink the size of an array an arbitrary amount
;   and compute uncertainty map for shrunk array
;   Based on idl astrolib FREBIN
;   Does not conserve flux, returns only mean maps.
;   ot tested on 1D arrays
;
;
; CALLING SEQUENCE:
;   SREBIN, image, dmap, nsout, nlout
;
; INPUTS:
;    image - input image, 1-d or 2-d numeric array
;    nsout - number of samples in the output image, numeric scalar
;
; OPTIONAL INPUT:
;    nlout - number of lines in the output image, numeric scalar
;            If not supplied, then set equal to 1
;
; OUTPUTS:
;    The resized image is returned as the function result.    If the input
;    image is of type DOUBLE or FLOAT then the resized image is of the same
;    type.     If the input image is BYTE, INTEGER or LONG then the output
;    image is usually of type FLOAT. 
;    The uncertainty map provides the error of the mean for each pixel.
;     
; EXAMPLE:
;     Suppose one has an 800 x 800 image array, im, that must be expanded to
;     a size 850 x 900 while conserving surface flux:
;
;     IDL> im1 = frebin(im,850,900,/total) 
;
;     im1 will be a 850 x 900 array, and total(im1) = total(im)
; 
; NOTES:
;  Rebinning factor MUST be >= 2, smaller values will probably lead to
;  error.
; 
; PROCEDURE CALLS:
;    None.
; HISTORY:
;    Adapted from May 1998 STIS  version, written D. Lindler, ACC
;    Added /NOZERO, use INTERPOLATE instead of CONGRID, June 98 W. Landsman  
;    Fixed for nsout non-integral but a multiple of image size  Aug 98 D.Lindler
;    DJL, Oct 20, 1998, Modified to work for floating point image sizes when
;		expanding the image. 
;    Improve speed by addressing arrays in memory order W.Landsman Dec/Jan 2001
;-
;    Modified by R. Schoedel to create uncertainty maps
;----------------------------------------------------------------------------

; determine size of input image
;
	ns = n_elements(image[*,0])
	nl = n_elements(image)/ns
	ns1 = ns-1
	nl1 = nl-1

        dtype = size(image,/TNAME)
	if dtype EQ 'DOUBLE' then begin
		sbox = ns/double(nsout) 
		lbox = nl/double(nlout)
	   end else begin
		sbox = ns/float(nsout) 
		lbox = nl/float(nlout)
	end	

; Now do 2-d case
; First, bin in second dimension
;
	    if dtype eq 'DOUBLE' then temp = dblarr(ns,nlout) $
			         else temp = fltarr(ns,nlout)

; loop on output image lines
;
	    for i=0L,nlout-1 do begin
	    	    rstart = i*lbox		;starting position for each box
	    	    istart = long(rstart)
	    	    rstop = rstart + lbox	;ending position for each box
	    	    istop = long(rstop)<nl1
	    	    frac1 = rstart-istart
	    	    frac2 = 1.0 - (rstop-istop)
;
; add pixel values from istart to istop and  subtract fraction pixel 
; from istart to rstart and fraction pixel from rstop to istop
;

                     if istart EQ istop then $
	   	       temp[0,i] = (1.0 - frac1 - frac2)*image[*,istart] $
                       else $
	   	       temp[0,i] = total(image[*,istart:istop],2) $
	   			- frac1 * image[*,istart]  $
	   			- frac2 * image[*,istop] 
	    endfor
            temp = transpose(temp)
;
; bin in first dimension
;
	    if dtype eq 'DOUBLE' then result = dblarr(nlout,nsout) $
			         else result = fltarr(nlout,nsout)

;
; loop on output image samples
;
	    for i=0L,nsout-1 do begin
	    	    rstart = i*sbox	       ;starting position for each box
	    	    istart = long(rstart)
	    	    rstop = rstart + sbox      ;ending position for each box
	    	    istop = long(rstop)<ns1
	    	    frac1 = rstart-istart
	    	    frac2 = 1.0 - (rstop-istop)
;
; add pixel values from istart to istop and  subtract fraction pixel 
; from istart to rstart and fraction pixel from rstop to istop
;

		    if istart eq istop then $
                        result[0,i] = (1.-frac1-frac2)*temp[*,istart] else $
		    	result[0,i] = total(temp[*,istart:istop],2)   $
		    		- frac1 * temp[*,istart]  $
		    		- frac2 * temp[*,istop]
	    end

;            
; return rebinned map
	    rmap = transpose(result)/(sbox*lbox)
	   
; now compute uncertainty map
; ==================================

; create mean map of same size as input image
meanim = FREBIN(rmap,ns,nl)

; Now do 2-d case
; First, bin in second dimension
;
	    if dtype eq 'DOUBLE' then temp = dblarr(ns,nlout) $
			         else temp = fltarr(ns,nlout)

; loop on output image lines
;
	    for i=0L,nlout-1 do begin
	    	    rstart = i*lbox		;starting position for each box
	    	    istart = long(rstart)
	    	    rstop = rstart + lbox	;ending position for each box
	    	    istop = long(rstop)<nl1
	    	    frac1 = rstart-istart
	    	    frac2 = 1.0 - (rstop-istop)


                     if istart EQ istop then $
	   	       temp[0,i] = (1.0 - frac1 - frac2)*(image[*,istart]-meanim[*,istart])^2 $
                       else $
	   	       temp[0,i] = total((image[*,istart:istop]-meanim[*,istart:istop])^2,2) $
	   			- frac1 * (image[*,istart]-meanim[*,istart])^2  $
	   			- frac2 * (image[*,istop] -meanim[*,istop])^2
	    endfor
            temp = transpose(temp)
;
; bin in first dimension
;
	    if dtype eq 'DOUBLE' then result = dblarr(nlout,nsout) $
			         else result = fltarr(nlout,nsout)

;
; loop on output image samples
;
	    for i=0L,nsout-1 do begin
	    	    rstart = i*sbox	       ;starting position for each box
	    	    istart = long(rstart)
	    	    rstop = rstart + sbox      ;ending position for each box
	    	    istop = long(rstop)<ns1
	    	    frac1 = rstart-istart
	    	    frac2 = 1.0 - (rstop-istop)
;
; add pixel values from istart to istop and  subtract fraction pixel 
; from istart to rstart and fraction pixel from rstop to istop
;

		    if istart eq istop then $
                        result[0,i] = (1.-frac1-frac2)*temp[*,istart] else $
		    	result[0,i] = total(temp[*,istart:istop],2)   $
		    		- frac1 * temp[*,istart]  $
		    		- frac2 * temp[*,istop]
	    end

;            
; return rebinned map
	    dmap = sqrt(transpose(result)/(sbox*lbox-1))  ; sigma
            dmap = dmap/sqrt(sbox*lbox)                   ; sigma --> uncertainty of the mean

end
