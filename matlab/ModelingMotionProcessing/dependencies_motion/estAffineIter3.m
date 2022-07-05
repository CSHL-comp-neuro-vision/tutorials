function M = estAffineIter3(vol1,vol2,numIters,Minitial)
%
% function M = estAffineIter3(vol1,vol2,numIters,Minitial)
%
% vol1 and vol2 are volumes, 3d arrays
% numIters is number of iterations to run
%
% M is 4x4 affine transform matrix: X' = M X
% where X=(x,y,z,1) is starting position in homogeneous coords
% and X'=(x',y',z',1) is ending position
%
% Each iteration warps the volumes according to the previous
% estimate, and estimates the residual motion.

if ~exist('numIters')
  numIters=3;
end

if ~exist('Minitial')
  %M=estAffine3(vol1,vol2);
  M=eye(4);
else
  M=Minitial;
end
%disp(['iter=0']);
%disp(M);

for iter=1:numIters
  Mhalf2=real(sqrtm(M));
  Mhalf1=real(sqrtm(inv(M)));
  volWarp1=warpAffine3(vol1,Mhalf1);
  volWarp2=warpAffine3(vol2,Mhalf2);
  disp(['iter=',num2str(iter)]);
  deltaM=estAffine3(volWarp1,volWarp2);
  M=deltaM*M;
  disp(M);
end

return;

%%%%%%%%%
% Debug %
%%%%%%%%%

% input
filter = [0.03504 0.24878 0.43234 0.24878 0.03504];
in = convXYsep(convZ(rand(30,40,14),filter),filter,filter);

% translation
A= [1 0 0 .2;
    0 1 0 .3;
    0 0 1 .4;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
vol1=vol1(:,:,2:9);
vol2=vol2(:,:,2:9);
A1*A1
Aest=estAffineIter3(vol1,vol2,3,eye(4))

% rotation in x,y
theta=.03;
A1=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0
    0 0 0 1];
A2=inv(A1);
vol1=warpAffine3(in,A1);
vol2=warpAffine3(in,A2);
vol1=vol1(:,:,2:9);
vol2=vol2(:,:,2:9);
A1*A1
M=estAffineIter3(vol1,vol2,5,eye(4));

% rotation in x,z
theta=.03;
A1=[cos(theta) 0 sin(theta) 0;
    0 1 0 0;
    -sin(theta) 0 cos(theta) 0;
    0 0 0 1];
A2=inv(A1);
vol1=warpAffine3(in,A1);
vol2=warpAffine3(in,A2);
vol1=vol1(:,:,2:9);
vol2=vol2(:,:,2:9);
A1*A1
M=estAffineIter3(vol1,vol2,5,eye(4));
