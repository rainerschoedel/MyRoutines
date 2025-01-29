PRO SORT_CUBES, maskdir, masklist, cube_order

  n_cubes = n_elements(masklist)

  mcube = readfits(maskdir + masklist[0])
  sz = size(mcube)
  n1 = sz[1]
  n2 = sz[2]
  fovs_sizes = lonarr(n_cubes)
  fovs = lonarr(n1,n2,n_cubes)
;  masks = lonarr(n1,n2,n_cubes)

  for ic = 0, n_cubes-1 do begin
    mcube = readfits(maskdir + masklist[ic])
    fov = mcube[*,*,0]  
;    writefits, 'tmp/fov.fits', fov
    ; We have to fill the holes in the fov,
    ; First, find the large parts of the fov
    ; near the edges that
    ; may zero. This will speed up thngs a lot.
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
    check = where(tmpfov lt 1,counts)
    for i = 0L, counts-1 do begin
      xy = array_indices(tmpfov,check[i])
      xx = xy[0]
      yy = xy[1]
      region = search2d(tmpfov,xx,yy,0,0.1)
      n_zero = n_elements(region)
      if (n_zero gt 0) then fov[region] = 1
    endfor
;    writefits, 'tmp/fov_filled.fits', fov
    valid = where(fov gt 0, fov_size)
    fovs_sizes[ic] = fov_size
    fovs[*,*,ic] = fov
  endfor
  cube_order = REVERSE(SORT(fovs_sizes))
  masklist = masklist[cube_order]
  print, masklist
  fovs = fovs[*,*,cube_order]
;  writefits, 'tmp/fovs_sorted.fits', fovs
;STOP  
  fov = replicate(1,n1,n2)
  for ic = 0, n_cubes-1 do begin
    mcube = readfits(maskdir + masklist[ic])
    thismask = mcube[*,*,0]  
;    writefits, 'tmp/thismask.fits', thismask
;    writefits, 'tmp/fov.fits', fov
    fov = fovs[*,*,ic] * fov
    fovs[*,*,ic] = fov
;    writefits, 'tmp/fov_masked.fits', fov
;    writefits, 'tmp/newmask.fits', fov * thismask
    sz = size(mcube)
    n3 = sz[3]
    thismask = thismask * fov
;    masks[*,*,ic] = thismask
    for i_m = 0, n3-1 do mcube[*,*,i_m] = thismask
    outname = 'holo' + STRMID(masklist[ic],0,STRLEN(masklist[ic])-3)
    masklist[ic] = 'holo' + masklist[ic]
    writefits, maskdir + outname, mcube, /COMPRESS
;    STOP
  endfor
;  writefits, 'tmp/fovs_masked.fits', fovs
;  writefits, 'tmp/masks.fits', masks
;  STOP

END
