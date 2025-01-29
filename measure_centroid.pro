PRO MEASURE_CENTROID, im, x, y, r0

 ; PURPOSE: Given a list x and y of positions in an image and a radius,
 ; compute the centroids at these positions within a given radius

 ; x and y must be FLOAT!

 xint = long(x)
 yint = long(y) 

 npos = n_elements(x)
 boxsize = 2*r0 + 1

 ; extract sub-arrays
 ; centers of sub-arrays will be [r0,r0]
 SUB_ARRAYS, im, xint, yint, boxsize, stack, masks

 for i = 0, (npos-1) do begin
  sub = stack[*,*,i]
  sub = circ_mask(sub,r0,r0, r0)
  cen = centroid(sub)
  x[i] = xint[i] + (cen[0]-r0)
  y[i] = yint[i] + (cen[1]-r0)
;  test 
;  print, x[i], y[i]
;  hhh = circ_mask(im,xint[i],yint[i],r0)
;  print, centroid(hhh)
;  STOP
 endfor

 

END
