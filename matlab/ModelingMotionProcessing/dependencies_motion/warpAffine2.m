function result = warpAffine2(im,A)
%
% function result = warpAffine3(im,A)
%
% im: input image
% A: 2x3 affine transform matrix or a 3x3 matrix with [0 0 1]
% for the last row.
% if a transformed point is outside of the volume, NaN is used
%
% result: output image, same size as im
%
% Author: David Heeger
%
% 8/13/97  dhb  Deleted extra 0 in comment above.

if (size(A,1)>2)
  A=A(1:2,:);
end

% Compute coordinates corresponding to input 
% and transformed coordinates for result
[x,y]=meshgrid(1:size(im,2),1:size(im,1));
coords=[x(:)'; y(:)'];
homogeneousCoords=[coords; ones(1,prod(size(im)))];
warpedCoords=A*homogeneousCoords;
xprime=warpedCoords(1,:)';
yprime=warpedCoords(2,:)';

result = interp2(x,y,im,xprime,yprime);
result = reshape(result,size(im));

return;

%%% Debug

im=[1 2 3; 4 5 6; 7 8 9]';

A= [1 0 .5;
    0 1 0;
    0 0 1];

A= [1 0 0;
    0 1 .5;
    0 0 1];

A= [1 0 .5;
    0 1 .5;
    0 0 1];

res=warpAffine2(im,A)


