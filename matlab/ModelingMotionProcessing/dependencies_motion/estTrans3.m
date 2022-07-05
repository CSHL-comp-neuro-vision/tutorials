function theta = estTrans3(in1,in2,dims)
%
% function theta = estTrans3(in1,in2,[m n])
%
% in1 and in2 are volumes, each column is an image
% [m n] is image size
%
% theta is 3-vector, translation in each (x,y,z) direction

[fx,fy,fz,ft]=computeDerivatives(in1,in2,dims);

A = [fx(:),fy(:),fz(:)];
b = -ft(:);

theta=A\b;

return;


% DEBUG:
dims=[20 20];
in1=rand(prod(dims),9);
%theta=[1 0 0];
%theta=[0 1 0];
theta=[0 0 1];
in2=circularShift3(in1,dims,theta);
estTrans3(in1,in2,dims)

