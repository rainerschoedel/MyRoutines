PRO LXP, indir, inlist, outim

; read names of input cubes
; -------------------------
name = ''
openr, inp, inlist, /get_lun
ncubes = 0L
while (not EOF(inp)) do begin
  readf, inp, name
  ncubes = ncubes + 1
endwhile
list = strarr(ncubes)
point_lun, inp, 0
for ic = 0, ncubes -1 do begin
  readf, inp, name
  list[ic] = name
endfor
free_lun, inp
print
print, 'input file list: ' 
print, list
print

; determine size of axes
; ------------------------
cube = readfits(indir + list[0])
sz = size(cube)
nax1 = sz[1]
nax2 = sz[2]

; initialize variables
; ---------------------

lxp = fltarr(nax1,nax2)
nt = 0L

for ic = 0, ncubes-1 do begin  ; start loop over all input cubes

; initialize variables
; ---------------------
  this = fltarr(nax1,nax2)
  nthis = 0L

  ; read speckle data cube
  ; --------------------------
  cube = readfits(indir + list[ic])
  cube = cube
  sz = size(cube)
  if (sz[0] gt 2) then  nim = sz[3] else nim = 1
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  print, 'Read in cube: ' + list[ic]
  print, 'Number of frames in this cube: ' + strtrim(string(nim),2)
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'

  for i = 0L, nim-1 do begin  ; start loop over frames in this cube

    im = cube[*,*,i]
    this = this + im
    nthis++
 
 endfor  ; end loop over frames in this cube

  writefits, '../tmp/lxp_' + strtrim(string(ic+1),2)+'.fits', this/nthis
  lxp = lxp + this
  nt = nt + nthis

endfor ; end loop over input cubes

print, 'Number of frames used: ' + strtrim(string(nt),2)
writefits, outim, lxp/nt

END
