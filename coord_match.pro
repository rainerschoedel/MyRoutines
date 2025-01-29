FUNCTION lintrans, P, X=X, Y=Y, ERR=err

 sz=size(X,/dimensions)
 ymod = dblarr(2,sz[1])
 ymod[0,*]=P[0]+P[2]*cos(P[3])*X[0,*]-P[2]*sin(P[3])*X[1,*]
 ymod[1,*]=P[1]+P[2]*sin(P[3])*X[0,*]+P[2]*cos(P[3])*X[1,*]
 
 D=sqrt((ymod[0,*]-y[0,*])^2+(ymod[1,*]-y[1,*])^2)/err
 D=transpose(D)
 return, D

END


PRO ITERPROC, lintrans, p, iter, fnorm, FUNCTARGS=fcnargs,_EXTRA=_extra

 XX=fcnargs.x
 YY=fcnargs.y
 error=fcnargs.err
 print, 'iteration #',iter
 DD=lintrans(p,X=XX,Y=YY,ERR=error)
 print, total(DD^2)
 print, p
END         

pro coord_match, x1, y1, x2, y2, dx1, dy1, params, Perror

;procedure for matching two lists of coordinates
;it allows for translation, rotation, and global difference in plate scales

;the form of the transformation is the following:
;
; x2=P[0]+P[2]*cos(P[3])*x1-P[2]*sin(P[3])*y1
; y2=P[1]+P[2]*sin(P[3])*x1+P[2]*cos(P[3])*y1
;
; i.e., [x1,y1] are rotated clockwise by P[3]
;
;where P[0]= offset along x
;      P[1]= offset along y
;      P[2]= ratio of the plate scales
;      P[3]= rotation angle

; x1, y1 --> X: positions in current frame that should be derotated
; x2, y2 --> Y: reference position(ssa image)
; dx1, dy1 --> Y: uncertainties in current frame


X=dblarr(2,n_elements(x1))
X[0,*]=x1
X[1,*]=y1
Y=dblarr(2,n_elements(x2))
Y[0,*]=x2
Y[1,*]=y2
err=dblarr(n_elements(x1))
err[*] = sqrt(dx1^2 + dy1^2)

parinfo = replicate({value:0.D, fixed:0,tied:''},4)

;starting values
parinfo[0].value=0.0
parinfo[1].value=0.0
parinfo[2].value=1.0
parinfo[3].value=0.0
parinfo[2].fixed=1  ; fix relative pixel scale

functargs = {X:X, Y:Y, ERR:err}

params = MPFIT('lintrans', FUNCTARGS=functargs, PARINFO=parinfo, PERROR=Perror, $
                          BESTNORM=chi2,ERRMSG=errmsg, MAXITER=1000)
;parms = MPFIT('lintrans', FUNCTARGS=functargs, PARINFO=parinfo, PERROR=Perror, $
;                          BESTNORM=chi2,ERRMSG=errmsg,ITERPROC='iterproc', MAXITER=1000)

DOF = n_elements(x1) - n_elements(params) ; deg of freedom
Perror = Perror * sqrt(chi2/DOF)   ; scaled uncertainties

end


