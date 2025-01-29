; $Id: repair_saturated.pro, v 1.1 Jul 2000 e.d. $
;
;
; NAME:
;	MYREPAIR_SATURATED
;
; PURPOSE:
;	Repair saturated stars in an image, by replacing their core with a
;	template representing an estimate of the PSF. Accurate positioning
;	of the template onto each star is performed by correlatio
;	maximization, whereas the scaling factor is computed by fitting the
;	wings of the saturated source.
;	This procedure gives acceptable results when the size of the
;	saturated region of a given star does not exceed the diameter of
;	the central core of the PSF.
;       This is a modified version of the original StarFinder code,
;       which does not work well when one has a limited knowledge of
;       the PSF wings.
;       The main point is to include an additive offset to the PSF.
;
; CATEGORY:
;	Signal processing. Stellar photometry.
;
; CALLING SEQUENCE:
;	REPAIR_SATURATED, Image, Clean_image, Background, $
;	                  Psf, Psf_fwhm, X, Y, Upper_lev
;
; INPUTS:
;	Image:	Stellar field image containing saturated stars.
;
;	Clean_image:	Stellar field image, after subtraction of non-saturated
;		sources. It may coincide with Image, if there are no important
;		sources other than the saturated ones.
;
;	Background:	2D array, with the same size as Image, containing an
;		estimate of the background emission.
;
;	Psf:	2D array, containing an estimate of the PSF, to be used as a
;		template to repair the core of the saturated stars.
;
;	Psf_fwhm:	FWHM of the PSF.
;
;	X, Y:	X- and Y- coordinates of saturated stars.
;
;	Upper_lev:	Scalar, representing the presumed saturation threshold.
;
; KEYWORD PARAMETERS:
;
;	N_WIDTH: Size of the box extracted around a saturated star is
;                N_WIDTH * widths[star]
;
;	MAG_FAC:	Integer representing the fractional sub-pixel step for
;		accurate positioning of the PSF estimate on the core of a saturated
;		star. This "magnification factor" is also used to optimize the
;		correlation.
;		The default is MAG_FAC = 2, corresponding to a positioning accuracy
;		of 1/2 pixel.
;
;	INTERP_TYPE:	Use this keyword to choose an interpolation technique
;		for the PSF fractional shift when MAG_FAC > 1. See IMAGE_SHIFT in
;		the file "image_shift.pro" for more details. The default is to
;		use cubic convolution interpolation.
;
; OUTPUTS:
;	Image:	Image array with repaired saturated stars.
;
;	Clean_image:	Input Clean_image with repaired saturated stars.
;
;	X, Y:	Coordinates of saturated stars after repair.
;
; SIDE EFFECTS:
;	The input variables Image, Clean_image, X and Y are overwritten.
;
; RESTRICTIONS:
;	It is assumed that the saturated stars are separated by a distance
;	at least (N_WIDTH * Width), where N_WIDTH is the input keyword described
;	above and Width is the maximum diameter in pixels of a saturated core.
;
; PROCEDURE:
;	The input image is temporarily cleaned from the contamination of
;	secondary sources and background emission. Then the saturated sources
;	are isolated and repaired with a scaled replica of the PSF template.
;	Accurate positioning of the template onto each star is performed by
;	correlation maximization, whereas the scaling factor is computed by
;	fitting the wings of the saturated source. Of course the saturated core
;	is masked when computing the correlation coefficient and the scaling
;	factor.
;	When all the saturated sources specified on input have been repaired,
;	the image is restored adding the previously subtracted stars and
;	the background emission.
;
; MODIFICATION HISTORY:
; 	Written by:	Emiliano Diolaiti, August-September 1999.
;	Updates:
;	1) Fixed bug on correlation maximization in module MATCH_REPLACE
;	   (Emiliano Diolaiti, July 2000).
;       2) Included additive offset
;          (Rainer Schoedel, March 2017)
;-




;;; Auxiliary routines.

PRO PSFFIT, X, P, YMOD
  
  YMOD = P[0] + P[1]*image_shift(X,P[2],P[3])

END


; MATCH_REPLACE: repair the core of a saturated star,
; given an unsaturated template.

PRO match_replace, image, template, x, y, box, $
                   sat_mask, MAG_FAC = mag_fac, x_fit = x_fit, $
                   y_fit = y_fit, f_fit = f_fit, _EXTRA = extra

        on_error, 2
        star = sub_array(image, box, REF = [x, y], $
					 LX = lx, UX = ux, LY = ly, UY = uy)
       ; Are there any NaNs present? --> set to 0 in W and star
        isnan = where(FINITE(star, /NAN))
        star[isnan] = 0
        writefits, 'star.fits', star
	sat = sat_mask[lx:ux,ly:uy]
        temp = sub_array(template, box)
        w = where(sat lt 1)
        P = [0.,1.,0.,0.] 
        P[1] = total(star)
        WT = star
        WT[*,*] = 0
        WT[w] = 1
        WT[isnan] = 0
        parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 4)
        parinfo[0].value = P[0]
        parinfo[1].value = P[1]
        parinfo[2].value = P[2]
        parinfo[3].value = P[3]
        parinfo[0].limited = [1,0]
        parinfo[1].limited = [1,0]
        parinfo[2].limited = [1,1]
        parinfo[3].limited = [1,1]
        parinfo[0].limits[0] = 0.
        parinfo[1].limits[0] = 0.
        parinfo[2].limits = [-box/4.,box/4.]
        parinfo[3].limits = [-box/4.,box/4.]

        res = mpcurvefit(temp,star,WT,P,sigma,FUNCTION_NAME='PSFFIT',/NODERIVATIVE,PARINFO=parinfo) ;,/QUIET)
        print, P

        ; Repair image
;        writefits, 'sat.fits', long(sat) ; illegal data type if I do not use long (???)
        writefits, 'W.fits', WT
        writefits, 'psf_template.fits', temp
        temp = P[0] + P[1] * image_shift(temp,P[2],P[3])
        writefits, 'psf_template_scaled.fits', temp
        writefits, 'diff.fits', WT * (star - temp)
	w = where(sat gt 0)
        WT = star
        WT[*,*] = 0
        WT[w] = 1
        writefits, 'toreplace.fits', WT
        star_sat = star
        star[w] = temp[w]
        writefits, 'star_rep.fits', star
        image[lx:ux,ly:uy] = star
        x_fit = x + P[2]
        y_fit = y + P[3]
        f_fit = P[1]
;STOP
        return
end



;;; The main routine.

PRO myrepair_saturated, image, clean_image, background, psf, x, y, widths, $
			sat_mask, N_WIDTH=n_width, sat_image = sat_image, _EXTRA = extra


        SetDefaultValue, n_width, 1
        on_error, 2

	; The image to use to repair saturated stars must be:
	; 1) background removed
        ; 2) cleaned from secondary stars around saturated ones.
        ; Careful with option 2: Do not clean
        ; PSF features around saturated star 
        sz = size(image) & n1 = sz[1] & n2 = sz[2]
        sec_sources = image - clean_image
        ; set NANs in sec_sources to 0 -->>
        ; those pixels will be masked later anyway
        isnan = where(FINITE(image,/NAN))
        sec_sources[isnan] = 0
        clean_image = temporary(clean_image) - background
	; Match and repair saturated stars in clean_image
        n_satur = n_elements(x)
        x_satur = fltarr(n_satur)
        y_satur = fltarr(n_satur)
        f_satur = fltarr(n_satur)
        n_fitted = 0L
	for  n = 0L, n_satur - 1  do begin
	   x_n = x[n]  &  y_n = y[n]
           print, 'Working on star '+ strn(n) + ' at ' + strn(x_n) + ', ' + strn(y_n) +'.'
           box = round(n_width * widths[n])
           box_hw = box/2
           ; Only repair a star if it is lies
           ; >box/2 from the image edges
           if (x_n gt (box_hw+1)) and (x_n lt (n1-box_hw-1)) and (y_n gt (box_hw+1)) and (y_n lt (n2-box_hw-1)) then begin
              match_replace, clean_image, psf, x_n, y_n, box, sat_mask, $
                             x_fit = x_fit, y_fit = y_fit, f_fit = f_fit, _EXTRA = extra
              print, 'Repaired star '+ strn(n) + ' at ' + strn(x_n) + ', ' + strn(y_n) +'.'
              x_satur[n_fitted] = x_fit
              y_satur[n_fitted] = y_fit
              f_satur[n_fitted] = f_fit
              n_fitted = n_fitted + 1
           endif else print, 'Star ' + strn(n) + ' is too close to edge: ' + strn(x_n) + ', ' + strn(y_n) + ', box: ' + strn(box) + '.'
;          STOP
        endfor
        x_satur = x_satur[0:n_fitted-1]
        y_satur = y_satur[0:n_fitted-1]
        f_satur = f_satur[0:n_fitted-1]
        
	; Define clean_image = input image with repaired stars - secondary sources
	clean_image = temporary(clean_image) + background
	; Define image = original image with repaired stars
        image = clean_image + sec_sources
;        writefits, 'clean_image.fits', clean_image

        ; make image of saturated stars only
         sat_image = image_model(x_satur,y_satur,f_satur,n1,n2,psf)
     
        return
end
