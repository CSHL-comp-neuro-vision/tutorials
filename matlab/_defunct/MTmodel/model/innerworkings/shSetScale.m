% res = shSetScale(shMatrix, ind, valuesToInsert, scale)    
%
% Replace all the values at one spatial scale in the matrix shMatrix with
% the entries in valuesToInsert.

function res = shSetScale(shMatrix, ind, valuesToInsert, scale)

scale = 1;
valuesToInsert = reshape(valuesToInsert, [size(valuesToInsert,1)*size(valuesToInsert,2)*size(valuesToInsert,3), size(valuesToInsert,4)]);
valuesToInsert = valuesToInsert';
shMatrix(:,ind(scale, 1)+1:ind(scale+1, 1)) = valuesToInsert;
res = shMatrix;
