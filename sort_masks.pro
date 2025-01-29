PRO SORT_MASKS, mcube, cube

  sz = size(mcube)
  n1 = sz[1]
  n2 = sz[2]
   n_masks = sz[3]

  fovs_sizes = lonarr(n_masks)
  fovs = lonarr(n1,n2,n_masks)

  for ic = 0, n_masks-1 do begin
    fov = mcube[*,*,ic]  
;    writefits, 'tmp/fov.fits', fov
    ; We have to fill the holes in the fov,
    ; First, find the large parts of the fov
    ; near the edges that
    ; may be zero. This will speed up things a lot.
    tmpfov = fov
    for x_edge = 0, 1 do begin
       for y_edge = 0, 1 do begin
         x0 = x_edge * n1
         y0 = y_edge * n2
         if (x_edge gt 0) then x0 = x0-1
         if (y_edge gt 0) then y0 = y0-1
         if (fov[x0,y0] eq 0) then begin
          region = search2D(fov,x0,y0,0,0.1)
          n_zero = n_elements(region)
          if (n_zero gt 0)then tmpfov[region] = 1
         endif
       endfor
    endfor
;    writefits, 'tmp/tmpfov.fits', tmpfov
    check = where(tmpfov lt 1,counts)
    if (counts gt 0) then fov[check] = 1
;    writefits, 'tmp/fov_filled.fits', fov
    valid = where(fov gt 0, fov_size)
    fovs_sizes[ic] = fov_size
    fovs[*,*,ic] = fov
    print, ic
;    STOP
  endfor
  mask_order = REVERSE(SORT(fovs_sizes))
  print, mask_order
;  fovs = fovs[*,*,mask_order]
;  writefits, 'tmp/fovs_sorted.fits', fovs
;STOP  
  fov = replicate(1,n1,n2)
  for ic = 0, n_masks-1 do begin
    thisfov = fovs[*,*,mask_order[ic]]
    thisim = cube[*,*,mask_order[ic]]  
    thismask = mcube[*,*,mask_order[ic]]  
;    writefits, 'tmp/thismask.fits', thismask
;    writefits, 'tmp/fov.fits', fov
     fov = thisfov * fov
;    writefits, 'tmp/fov_masked.fits', fov
;    writefits, 'tmp/newmask.fits', fov * thismask
    thismask = thismask * fov
    mcube[*,*,ic] = thismask
    cube[*,*,ic] = thismask * thisim
;    STOP
    print, ic
  endfor

END
