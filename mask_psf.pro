PRO MASK_PSF, psf, maskrad, PSF_MASKED=psf_masked, PSF_OFFSET=psf_offset, WINGS=wings

; PURPOSE: A circular mask is applied to a PSF and the offset of the
;          PSF is subtracted. The offset is determined from a one pixel
;          wide ring immediately outside the masking radius.

sz = size(psf)
n1 = sz[1]
n2 = sz[2]
psf_size = n1
mid = long(psf_size/2)

ringmask = replicate(1, psf_size, psf_size)
ringmask = CIRC_MASK(ringmask, mid, mid, maskrad, /INNER)
;ringmask = CIRC_MASK(ringmask, mid, mid, maskrad+1)
off_ind = where(ringmask gt 0)
psf_offset = median(psf[off_ind])
print, 'PSF offset determined by MASK_PSF: ' + strn(psf_offset)

psf_masked = psf - psf_offset
psf_masked = CIRC_MASK(psf_masked, mid, mid, maskrad)
wings = CIRC_MASK(psf, mid, mid, maskrad, /INNER)



END
