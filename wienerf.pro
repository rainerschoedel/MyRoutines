PRO PSDMODEL, X, P, YMOD

  sz = size(X)
  n1 = sz[1]
  n2 = sz[2]
  xm = n1/2 
  ym = n2/2 
  YMOD = fltarr(n1,n2)
  for i = 0, n1 - 1 do begin
     for j = 0, n2-1 do begin
      YMOD[i,j] = P[0] + exp(-((i-xm)^2+(j-ym)^2)/(2*P[1]^2)) * P[2]
     endfor
  endfor

END

PRO WIENERF, image, params = params

  sz = size(image)
  n1 = sz[1]
  n2 = sz[2]
  H = FFT(image)
  HABS = ABS(H)
  HABS = shift(HABS,n1/2,n2/2)


  writefits, 'psd.fits', HABS
STOP
 ; fit model PSD to data: Gaussian + constant
  xm = float(n1/2.)
  ym = float(n2/2.)
  Y = HABS


  if (keyword_set(params)) then P = params else P = [0.,float(n1/2),median(Y)]
   X = Y
   X[*,*] = 1.
   W = X
   psdmodel = mpcurvefit(X,Y,W,P,sigma,FUNCTION_NAME='PSDMODEL',/NODERIVATIVE);,/QUIET)
 ;  print, P
 ;  writefits, 'model.fits', model
 ;  writefits, 'resid.fits', HABS-model

  N = fltarr(n1,n2)
  N[*,*] = P[0]
  S = fltarr(n1,n2)
  for i = 0, n1 - 1 do begin
     for j = 0, n2-1 do begin
      S[i,j] = exp(-((i-xm)^2+(j-ym)^2)/(2*P[1]^2)) * P[2]
     endfor
  endfor
  W = S/(S+N)
  W = shift(W,-n1/2,-n2/2)
;  writefits, 'HABS_filtered.fits', HABS * W

  O = H * W
  image = REAL_PART(fft(O,/INVERSE))

  params = P
;  writefits, 'wiener.fits', image

END
