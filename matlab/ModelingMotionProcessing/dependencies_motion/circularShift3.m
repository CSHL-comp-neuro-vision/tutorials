function result = circularShift3(input,theta)
%
% result = circularShift3(input,theta)
%
% theta is a 3-vector of integers

result=zeros(size(input));
tmp=zeros(size(input));

% circular shift in X and Y
for sliceNum=1:size(input,3)
  tmp(:,:,sliceNum)=circularShift(input(:,:,sliceNum),theta(1),theta(2));
end

% circular shift in Z
zShift=theta(3);
lastSlice=size(input,3);
if (zShift>0)
  result(:,:,[1:zShift]) = tmp(:,:,[lastSlice-zShift+1:lastSlice]);
  result(:,:,[zShift+1:lastSlice]) = tmp(:,:,[1:lastSlice-zShift]);
else
  zShift = -zShift;
  result(:,:,[1:lastSlice-zShift]) = tmp(:,:,[zShift+1:lastSlice]);
  result(:,:,[lastSlice-zShift+1:lastSlice]) = tmp(:,:,[1:zShift]);
end

return;

% Debug:
input=ones(3,3,4);
for z=1:4
  input(:,:,z)=z*input(:,:,z);
end
circularShift3(input,[0 0 1])
circularShift3(input,[0 0 -1])

slice=[1 2 3; 4 5 6; 7 8 9];
input=ones(3,3,4);
for z=1:4
  input(:,:,z)=slice;
end
result=circularShift3(input,[1 0 0])
result=circularShift3(input,[-1 0 0])
result=circularShift3(input,[0 1 0])
result=circularShift3(input,[0 -1 0])
