FUNCTION READLIST, listname, nc

line = fltarr(nc)

openr, inp, listname, /get_lun
nl = 0L
while (not EOF(inp)) do begin
  readf, inp, line
  nl++
endwhile
list = fltarr(nc,nl)
point_lun, inp, 0
for i = 0, nl-1 do begin
  readf, inp, line
  list[*,i] = line
endfor
free_lun, inp

return, list

END
