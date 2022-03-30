PRO RETURNMARKED, xsize, ysize, x, y, f, XM = xm, YM = ym, FM = fm, BOXSIZE = boxsize, dmax = dmax, g_sigma = g_sigma, DISP_STRETCH = disp_stretch, DISP_LARGE = disp_large, index = index, DISPRAN = dispran

;PURPOSE: Receive list of image size, stellar positions and fluxes
;         Create artifical map, dysplay it, and return marked stars.

device, decomposed = 0

if not(KEYWORD_SET(boxsize)) then boxsize = 5.
if not(KEYWORD_SET(dmax)) then dmax = 1.
if not(KEYWORD_SET(g_sigma)) then g_sigma = 1.5
if not(KEYWORD_SET(disp_stretch)) then disp_stretch = 'linear'
if not(KEYWORD_SET(disp_large)) then disp_large = 0
if not (KEYWORD_SET(dispran)) then dispran = [1.0e-5,0.1]


dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: g_sigma, Sigma_y: g_sigma, Angle: 0.0})
im = image_model(x,y,f,xsize,ysize,'gaussian', dat)

disp_opt = default_display_opt(im)
disp_opt.stretch = disp_stretch
disp_opt.large = disp_large
disp_opt.range = max(im)* dispran
display_image, im, wnum, OPTIONS = disp_opt, MODIFY_OPT = modify_opt

click_on_max, im, /MARK, BOXSIZE = boxsize, x_mark, y_mark

compare_lists, x_mark, y_mark, x, y, x1c, y1c, x2c ,y2c, SUBSCRIPTS_1 = subc1, SUBSCRIPTS_2 = subc2, MAX_DISTANCE=dmax
nc = n_elements(x1c)
print, 'Number of marked stars found in list: ' + string(nc)
xm = x[subc2]
ym = y[subc2]
fm = f[subc2]
index = subc2

WDELETE, wnum

END
