function result = warpAffine3(in,A,badVal)
%
% function result = warpAffine3(in,A,badVal)
%
% in: input volume, 3D array
% A: 3x4 affine transform matrix or a 4x4 matrix with [0 0 0 1]
%    for the last row.
% badVal: if a transformed point is outside of the volume, badVal is used
%
% result: output volume, same size as in

if ~exist('badVal')
  badVal=NaN;
end

if (size(A,1)>3)
  A=A(1:3,:);
end

% Compute coordinates corresponding to input volume
% and transformed coordinates for result
[xgrid,ygrid,zgrid]=meshgrid(1:size(in,2),1:size(in,1),1:size(in,3));
coords=[xgrid(:)'; ygrid(:)'; zgrid(:)'];
homogeneousCoords=[coords; ones(1,size(coords,2))];
warpedCoords=A*homogeneousCoords;

% Compute result using interp3
result = interp3(in,warpedCoords(1,:),warpedCoords(2,:),warpedCoords(3,:),'*linear*');
result = reshape(result,size(in));

% replace NaNs with badval
if(~isnan(badVal)) 
  NaNIndices = find(isnan(result));
  result(NaNIndices)=badVal*ones(size(NaNIndices));
end
return;

%%% Debug

slice=[1 2 3; 4 5 6; 7 8 9]';
slice=[1 1 1; 3 3 3; 5 5 5]';
input=ones(3,3,4);
for z=1:4
  input(:,:,z)=slice;
end

A= [1 0 0 .5;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

A= [1 0 0 .5;
    0 1 0 .5;
    0 0 1 0];

res=warpAffine3(input,A)
res=warpAffine3(input,A,-1)

for z=1:4
  input(:,:,z)=z*ones(3,3);
end
A= [1 0 0 0;
    0 1 0 0;
    0 0 1 .5;
    0 0 0 1];
res=warpAffine3(input,A)

input=rand(5,5,5);
res=warpAffine3(input,eye(4));
