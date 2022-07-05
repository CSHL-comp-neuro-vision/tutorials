function M = estAffine3(vol1,vol2)
%
% function M = estAffine3(vol1,vol2)
%
% vol1 and vol2 are volumes, 3d arrays
%
% M is 4x4 affine transform matrix: X' = M X
% where X=(x,y,z,1) is starting position in homogeneous coords
% and X'=(x',y',z',1) is ending position
%
% Solves fs^t theta + ft = 0
% where theta = B p is image velocity at each pixel
%       B is 3x12 matrix that depends on image positions
%       p is vector of affine motion parameters
%       fs is vector of spatial derivatives at each pixel
%       ft is temporal derivative at each pixel
% Mulitplying fs^t B gives a 1x12 vector for each pixel.  Piling
% these on top of one another gives A, an Nx12 matrix, where N is
% the number of pixels.  Solve M p = ftVol where ftVol is an Nx1
% vector of the the temporal derivatives at every pixel.

[fxVol,fyVol,fzVol,ftVol]=computeDerivatives3(vol1,vol2);
numSlices=size(fxVol,3);
[xgrid,ygrid,zgrid]=meshgrid(1:size(vol1,2),1:size(vol1,1),1:size(vol1,3));

% subsample
dims=size(fxVol);
fxVol=fxVol([1:2:dims(1)],[1:2:dims(2)],:);
fyVol=fyVol([1:2:dims(1)],[1:2:dims(2)],:);
fzVol=fzVol([1:2:dims(1)],[1:2:dims(2)],:);
ftVol=ftVol([1:2:dims(1)],[1:2:dims(2)],:);
% *** Assumes that the filters have 5 taps!
xgrid=xgrid([3:2:dims(1)+2],[3:2:dims(2)+2],[3:dims(3)+2]);
ygrid=ygrid([3:2:dims(1)+2],[3:2:dims(2)+2],[3:dims(3)+2]);
zgrid=zgrid([3:2:dims(1)+2],[3:2:dims(2)+2],[3:dims(3)+2]);

pts=find(~isnan(fxVol));
disp(['numPts=',num2str(length(pts))]);
fxVol = fxVol(pts);
fyVol = fyVol(pts);
fzVol = fzVol(pts);
ftVol = ftVol(pts);
xVol=xgrid(pts);
yVol=ygrid(pts);
zVol=zgrid(pts);

A= [xVol(:).*fxVol(:), yVol(:).*fxVol(:), zVol(:).*fxVol(:), fxVol(:),...
    xVol(:).*fyVol(:), yVol(:).*fyVol(:), zVol(:).*fyVol(:), fyVol(:),...
    xVol(:).*fzVol(:), yVol(:).*fzVol(:), zVol(:).*fzVol(:), fzVol(:)];
b = -ftVol(:);
p = A\b;

M= [1+p(1) p(2) p(3) p(4);
    p(5) 1+p(6) p(7) p(8);
    p(9) p(10) 1+p(11) p(12);
    0 0 0 1];

return;

%%%%%%%%%
% Debug %
%%%%%%%%%

% test with translation
vol1=rand(16,16,9);
%theta=[1 0 0];
%theta=[0 1 0];
theta=[0 0 1];
vol2=circularShift3(vol1,theta);
estAffine3(vol1,vol2)

in=rand(30,32,9);
A= [1 0 0 .2;
    0 1 0 .3;
    0 0 1 .4;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
Aest=estAffine3(vol1,vol2)

% test with rotation
in=rand(30,32,9);
theta=.01;
A1=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0;
    0 0 0 1];
A2=inv(A1);
vol1=warpAffine3(in,A1);
vol2=warpAffine3(in,A2);
A1*A1
Aest=estAffine3(vol1,vol2)

