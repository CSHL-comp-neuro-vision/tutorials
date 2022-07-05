function M = estAffine2(im1,im2)
%
% function M = estAffine2(im1,im2)
%
% im1 and im2 are images
%
% M is 3x3 affine transform matrix: X' = M X
% where X=(x,y,1) is starting position in homogeneous coords
% and X'=(x',y',1) is ending position
%
% Solves fs^t theta + ft = 0
% where theta = B p is image velocity at each pixel
%       B is 2x6 matrix that depends on image positions
%       p is vector of affine motion parameters
%       fs is vector of spatial derivatives at each pixel
%       ft is temporal derivative at each pixel
% Mulitplying fs^t B gives a 1x6 vector for each pixel.  Piling
% these on top of one another gives A, an Nx6 matrix, where N is
% the number of pixels.  Solve M p = ft where ft is now an Nx1
% vector of the the temporal derivatives at every pixel.

[fx,fy,ft]=computeDerivatives2(im1,im2);
[xgrid,ygrid]=meshgrid(1:size(im1,2),1:size(im1,1));

% subsample
dims=size(fx);
fx=fx([1:2:dims(1)],[1:2:dims(2)]);
fy=fy([1:2:dims(1)],[1:2:dims(2)]);
ft=ft([1:2:dims(1)],[1:2:dims(2)]);
% *** Assumes that the filters have 5 taps!
xgrid=xgrid([3:2:dims(1)+2],[3:2:dims(2)+2]);
ygrid=ygrid([3:2:dims(1)+2],[3:2:dims(2)+2]);


pts=find(~isnan(fx));
fx = fx(pts);
fy = fy(pts);
ft = ft(pts);
xgrid = xgrid(pts);
ygrid = ygrid(pts);

A= [xgrid(:).*fx(:), ygrid(:).*fx(:), fx(:),...
    xgrid(:).*fy(:), ygrid(:).*fy(:), fy(:)];
b = -ft(:);
p = A\b;

M= [1+p(1) p(2) p(3);
    p(4) 1+p(5) p(6);
    0 0 1];

return;

%%%%%%%%%
% Debug %
%%%%%%%%%

% test with translation
dims=[64 64];
im1=rand(dims);
im2=circularShift(im1,1,0);
%im2=circularShift(im1,0,1);
estAffine2(im1,im2)

dims=[64 64];
im1=rand(dims);
A= [1 0 .5;
    0 1 .5;
    0 0 1];
im2=warpAffine2(im1,A);
Aest=estAffine2(im1,im2)

% test with rotation
dims=[64 64];
in=rand(dims);
theta=atan2(1,max(dims));
A1=[cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];
A2=inv(A1);
im1=warpAffine2(in,A1);
im2=warpAffine2(in,A2);
A1*A1
Aest=estAffine2(im1,im2)

