FUNCTION SHIFTNW, im, xs, ys

; PURPOSE: Like IDL SHIFT procedure, but no wrapping of image around
; edges
; INPUT: image must be 2-dimensional
;        shifts must be integer

sz = size(im)
n1 = sz[1]
n2 = sz[2]

sim = shift(im,xs,ys)
if (xs gt 0) then sim[0:xs-1,*] = 0 $
 else if (xs lt 0) then sim[n1+xs:n1-1,*] = 0
if (ys gt 0) then sim[*,0:ys-1] = 0 $
 else if (ys lt 0) then sim[*,n2+ys:n2-1] = 0

return, sim

END
