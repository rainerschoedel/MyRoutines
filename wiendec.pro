FUNCTION WIENDEC, im, psf, NSUP = nsup

if (n_elements(nsup) lt 1) then nsup = 0.01

sz = size(im)
n1 = sz[1]
n2 = sz[2]
cen1 = n1/2
cen2 = n2/2
center = [cen1,cen2]

sz = size(psf)
pn1 = sz[1]
pn2 = sz[2]

p = fltarr(n1,n2)
p[0:pn1-1,0:pn2-1] = psf
xs = n1/2 - pn1/2
ys = n2/2 - pn2/2
psf = image_shift(p,xs,ys)

H = FFT(psf)*(n1*n2)
HQ = CONJ(H)
HABS = ABS(H)
H2 = HABS^2
S = FFT(im)
SABS = ABS(S)
D = H2
mH2 = moment(H2)
D[*,*] = nsup*mH2[0]
hhh = REAL_PART(FFT((S * HQ*HABS/(H2*HABS+D)), /INVERSE))
rcim = shift(hhh, center)

return, rcim

END
