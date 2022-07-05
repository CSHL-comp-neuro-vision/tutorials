function M = estAffineIter2(im1,im2,numIters,Minitial)
%
% function M = estAffineIter2(im1,im2,numIters,Minitial)
%
% im1 and im2 are input images
% numIters is number of iterations to run
% Minitial is initial guess for M.  Default is 3x3 identity matrix.
%
% M is 3x3 affine transform matrix: X' = M X
% where X=(x,y,1) is starting position in homogeneous coords
% and X'=(x',y'',1) is ending position
%
% Each iteration warps the images according to the previous
% estimate, and estimates the residual motion.

if ~exist('numIters')
  numIters=3;
end

if ~exist('Minitial')
  M=estAffine2(im1,im2);
else
  M=Minitial;
end
disp(['iter=0']);
disp(M);

for iter=1:numIters
  Mhalf2=real(sqrtm(M));
  Mhalf2=Mhalf2(1:3,:);
  Mhalf1=real(sqrtm(inv(M)));
  Mhalf1=Mhalf1(1:3,:);
  imWarp1=warpAffine2(im1,Mhalf1);
  imWarp2=warpAffine2(im2,Mhalf2);
  deltaM=estAffine2(imWarp1,imWarp2);
  M=deltaM*M;
  disp(['iter=',num2str(iter)]);
  disp(M);
end

return;

%%%%%%%%%
% Debug %
%%%%%%%%%

dims=[256 96];
in=rand(dims);
theta=atan2(1,max(dims));
A1=[cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];
A2=inv(A1);
im1=warpAffine2(in,A1);
im2=warpAffine2(in,A2);
A1*A1
M=estAffineIter2(im1,im2,10);
