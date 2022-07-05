function M = estAffineMulti2(im1,im2,iters,Minitial)
%
% function M = estAffineMulti2(im1,im2,iters,Minitial)
%
% im1 and im2 are input images
% iters is vector of number of iterations to run at each
%   successive scale. default is: [3] that runs 3 iterations at
%   the base scale.
% Minitial is initial guess for M.  Default is 3x3 identity matrix.
% 
% This function calls itself recursively.
%
% M is 3x3 affine transform matrix (in homogeneous coordinates)
%
% Bugs and limitations:
% - should check that images are big enough for requested number
%   of levels. 
%
% Author: David Heeger
%
% 7/30/97  dhb, eah  Special case iters(1) == 0 to allow 
%                    specification of not doing finest scales.

if ~exist('iters')
  iters=[3];
end
if ~exist('Minitial')
  Minitial=eye(3);
end

if (length(iters)>1)
  % reduce images
  im1Small=reduce(im1);
  im2Small=reduce(im2);

  % reduce inital affine matrix
  M=Minitial;
  M(1:2,3)=M(1:2,3)/2;

  % estimate M for reduced images
  M=estAffineMulti2(im1Small,im2Small,iters(2:length(iters)),M);

  % expand estimated affine matrix
  M(1:2,3)=M(1:2,3)*2;
  Minitial=M;
end

% Iterate, warping and refining estimates
if (iters(1) > 0)
  M=estAffineIter2(im1,im2,iters(1),Minitial);
else
  M = Minitial;
end

return;

%%%%%%%%%
% Debug %
%%%%%%%%%

% Test translation
dims=[128 128];
in=rand(dims);
A1=[1 0 2;
    0 1 2;
    0 0 1];
A2=inv(A1);
im1=warpAffine2(in,A1);
im2=warpAffine2(in,A2);
A1*A1
% One-scale only fails
M=estAffineIter2(im1,im2,3);
% Multiscale wins
M=estAffineMulti2(im1,im2,[3,1,1]);

% Test rotation
dims=[256 128];
in=rand(dims);
theta=4*atan2(1,max(dims));
A1=[cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];
A2=inv(A1);
im1=warpAffine2(in,A1);
im2=warpAffine2(in,A2);
A1*A1
M=estAffineIter2(im1,im2,3)
M=estAffineMulti2(im1,im2,[3,2,2]);

% Test expansion
dims=[128 128];
in=rand(dims);
s=1.05
A1=[s 0 0;
    0 s 0;
    0 0 1];
A2=inv(A1);
im1=warpAffine2(in,A1);
im2=warpAffine2(in,A2);
A1*A1
M=estAffineIter2(im1,im2,3)
M=estAffineMulti2(im1,im2,[3,2,2]);

