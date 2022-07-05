function M = estAffineMulti3(vol1,vol2,iters,Minitial)
%
% function M = estAffineMulti3(vol1,vol2,iters,Minitial)
%
% vol1 and vol2 are volumes, 3d arrays
% iters is vector of number of iterations to run at each
%   successive scale. default is: [3] that runs 3 iterations at
%   the base scale.
% Minitial is initial guess for M.  default is 4x4 identity matrix.
%
% This function calls itself recursively.
%
% Bugs and limitations:
% - currently supports subsampling only in x and y.
% - should check that images are big enough for requested number
%   of levels. 

if ~exist('iters')
  iters=[3];
end
if ~exist('Minitial')
  Minitial=eye(4);
end
M=Minitial;

numSlices=size(vol1,3);

if (length(iters)>1)
  % reduce images
  newdims=size(reduce(vol1(:,:,1)));
  vol1small=zeros(newdims(1),newdims(2),numSlices);
  vol2small=zeros(newdims(1),newdims(2),numSlices);
  for z=1:numSlices
    vol1small(:,:,z)=reduce(vol1(:,:,z));
    vol2small(:,:,z)=reduce(vol2(:,:,z));
  end
  % reduce affine matrix
  M(1:2,4)=M(1:2,4)/2;
  % estimate M for reduced images
  M=estAffineMulti3(vol1small,vol2small,iters(2:length(iters)),M);
  % expand estimated affine matrix
  M(1:2,4)=M(1:2,4)*2;
end

% Iterate, warping and refining estimates
M=estAffineIter3(vol1,vol2,iters(1),M);

return;


%%%%%%%%%
% Debug %
%%%%%%%%%

% input
filter = [0.03504 0.24878 0.43234 0.24878 0.03504];
in = convXYsep(convZ(rand(68,88,14),filter),filter,filter);

% translation
A1=[1 0 0 1;
    0 1 0 1.5;
    0 0 1 0.3
    0 0 0 1];
A2=inv(A1);
vol1=warpAffine3(in,A1);
vol2=warpAffine3(in,A2);
vol1=vol1(:,:,2:9);
vol2=vol2(:,:,2:9);
M=estAffineIter3(vol1,vol2,3,eye(4));
M=estAffineMulti3(vol1,vol2,[3,3,3]);
A1*A1

% rotation in x,y
theta=.06;
A1=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0
    0 0 0 1];
A2=inv(A1);
vol1=warpAffine3(in,A1);
vol2=warpAffine3(in,A2);
vol1=vol1(:,:,2:9);
vol2=vol2(:,:,2:9);
M=estAffineIter3(vol1,vol2,3);
M=estAffineMulti3(vol1,vol2,[3,3,3]);
A1*A1

% rotation in x,z
theta=.06;
A1=[cos(theta) 0 sin(theta) 0;
    0 1 0 0;
    -sin(theta) 0 cos(theta) 0;
    0 0 0 1];
A2=inv(A1);
vol1=warpAffine3(in,A1);
vol2=warpAffine3(in,A2);
vol1=vol1(:,:,2:9);
vol2=vol2(:,:,2:9);
M=estAffineIter3(vol1,vol2,3);
M=estAffineMulti3(vol1,vol2,[3,3,3]);
A1*A1

% expansion
s=1.05;
A1=[s 0 0 0;
    0 s 0 0;
    0 0 1/s 0;
    0 0 0 1];
A2=inv(A1);
vol1=warpAffine3(in,A1);
vol2=warpAffine3(in,A2);
vol1=vol1(:,:,2:9);
vol2=vol2(:,:,2:9);
M=estAffineIter3(vol1,vol2,3);
M=estAffineMulti3(vol1,vol2,[3,3,3]);
A1*A1
