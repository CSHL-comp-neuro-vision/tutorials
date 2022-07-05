function [err,Xwhiterot,R] = kurtMatErr(p,Xwhite)

R = [cos(p.ang), sin(p.ang) ; -sin(p.ang) cos(p.ang)];
Xwhiterot = Xwhite*R;
err = -(kurtMat(Xwhiterot)-3).^2;