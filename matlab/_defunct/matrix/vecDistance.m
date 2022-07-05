function d = vecDistance(X,Y)
%
%  d = vecDistance(X,Y)
%
%AUTHOR:  Wandell
%DATE:  07.07.95 
%PURPOSE:  
%
%  Compute the Euclidean distance between the corresponding columns
% of the matrices X and Y.  The location vectors, (x,y,z), are contained
% in the rows of X and Y respecctively.
%
%   If X and Y are both matrices, then the routine
% computes the Euclidean distance between the corresponding tows
% of X and Y.  
%
%   If X is a vector, then the routine computes
% the Euclidean distance between the vector x and each row of Y.
%
%ARGUMENTS:
%  X,Y:  Matrices (or vector and matrix) whose columns contain the
%        coordinates of a set of points
%RETURNS:
%  d:    vector of distances
%
%HISTORY:
% 8.22.95 -- Changed row/col structure for inputs, BW

nX = size(X,1);
nY = size(Y,1);

if nX > 1
 if nX ~= nY
   error('euclideanDistance:  matrices must have the same number of columns')
 end

 delta = X - Y;

elseif nX == 1
 delta = Y - X(ones(nY,1),:);
end

d = sqrt(diag(delta*delta'));
