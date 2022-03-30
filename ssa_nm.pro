PRO SSA_NM, indir, outdir, inlist, outim,  maskrad,  refsource, rebfac, debug, select, med_filt, OUTPUT = output

; Like ssa.pro, but without masks
; that need to be shifted.
; --------------------------------

name = ''
openr, inp, inlist, /get_lun
ic = 0L

; determine size of image
; ------------------------
readf, inp, name
point_lun, inp, 0
cube = readfits(indir + name, NSLICE=1)
sz = size(cube)
nax1 = sz[1]
nax2 = sz[2]

; multiply with rebin factor
; --------------------------
maskrad = rebfac*maskrad
refsource = rebfac*refsource
nax1 = rebfac*nax1
nax2 = rebfac*nax2
cen1 = nax1/2
cen2 = nax2/2
center = [cen1,cen2]
; initialize variables
; ---------------------
ssa = fltarr(nax1,nax2)
nt = 0L
if (select eq 0) then $
  openw, out, '../data/stat.txt', /get_lun

while (not EOF(inp)) do begin
  readf, inp, name

  ; initialize variables
  ;---------------------

  this = fltarr(nax1,nax2)
  nthis = 0L

  ; read speckle data cube
  ; --------------------------
  cube = readfits(indir + name)
  sz = size(cube)
  if (sz[0] gt 2) then  nim = sz[3] else nim = 1
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  print, 'Read in cube: ' + name
  print, 'Number of frames in this cube: ' + strtrim(string(nim),2)
  print, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'


  for i = 0L, nim-1 do begin  ; start loop over frames in this cube

    ; load ith image 
    im = cube[*,*,i]

    ; optionally rebin image
    if (rebfac gt 1) then begin
;     writefits, '../tmp/origim.fits', im
      im = CONGRID(im,nax1,nax2,/CENTER,CUBIC = -0.5)
    endif

; PSF extraction
; ################

    x_psf = refsource[0]
    y_psf = refsource[1]
    nref = n_elements(x_psf)
    xint = long(x_psf)
    yint = long(y_psf)
    sub_arrays, im, xint, yint, nax1, stack, masks
    psf = stack[*,*,0]
    psf = circ_mask(psf,cen1,cen2,maskrad,BORDER=0)
    if (med_filt gt 0) then psf = median_filter(psf,med_filt)
    psfmax = max(psf, max_subscript)
    print, string(i+1) + ' ' + string(psfmax)
    if (select eq 0) then printf, out, psfmax
    mind = array_indices(psf,max_subscript)
    if (debug gt 0) then begin
      writefits, '../tmp/psf.fits', psf
    endif    

; shift to max for SSA and to centroid for holo-prep
; #######################################################

    if (psfmax gt select) then begin
      xs = round(cen1 -mind[0])
      ys = round(cen2 - mind[1])
      if (debug gt 0) then begin
       writefits, '../tmp/im.fits', im
      endif
      hhh = shiftnw(im,xs,ys)  
      this = this + hhh
     ;centroid determination may go nuts if there are strong negativities
     ; set negative values = 0 to avoid this bug
     ; of course, this assumes that values< 0 are not important...
      neg = where(psf lt 0)
      if (neg[0] gt -1) then psf[neg] = 0
      centr = centroid(psf)
      xs = round(cen1 - centr[0])
      ys = round(cen2 - centr[1])
      cube[*,*,nthis] = shiftnw(im,xs,ys)
      if (debug gt 0) then begin
       print, 'shift: ' + string(xs) + ', ' + string(ys)
       writefits, '../tmp/psf.fits', psf
       writefits, '../tmp/im.fits', cube[*,*,nthis]
       STOP
     endif    
     nthis++
    endif

  endfor  ; end loop over frames in this cube

  ic = ic + 1
  if (nthis gt 0) then begin
    writefits, '../tmp/ssa_'+ strtrim(string(ic),2) + '.fits', this/nthis
    ssa = ssa + this
    nt = nt + nthis 
    if (output gt 0) then $
     writefits, outdir + name, cube[*,*,0:nthis-1]
  endif
  ; write selected frames to file

endwhile ; end loop over input cubes

free_lun, inp
if (select eq 0) then free_lun, out

print, 'Number of frames used: ' + strtrim(string(nt),2)
writefits, outim, ssa/nt

END

