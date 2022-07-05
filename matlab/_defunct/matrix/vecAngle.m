function [ theta ] = vecAngle(X,Y)
%
% [theta]=vecAngle(X,Y)
%
% AUTHOR:  B. Wandell
% PURPOSE:
%
%   If X and Y are both matrices then this routine computes
%   the angles between the corresponding row vectors in X and Y.
%
%   vecAngle(x,Y)
%   if x is a vector, then this routine computes the angle between
%   x and each of the rows of Y
%
%   The angle, theta, between two vectors is related to the dot product by
%
%      x . y =  ||x|| * ||y|| cos(theta)
%
%
% HISTORY:
% 8.22.95 -- Changed the row/col structure, BW

%DEBUG:
% X = rand(5,3);
% Y = rand(5,3);
% X = rand(1,3);
% X = eye(3,3);
% Y = eye(3,3)
% X = [1 0 0]

nX = size(X,1);
nY = size(Y,1);

if nX > 1

 if nX ~= nY
   error('vecAngle:  matrices must have the same number of columns')
 end

 innerProd = diag(X*Y');
 l = sqrt(diag(X*X')) .* sqrt(diag(Y*Y'));
 theta = acos(innerProd ./ l);
 
elseif nX == 1
 
 innerProd = X*Y';
 l = (norm(X,'fro'))*sqrt(diag(Y*Y'));

%The real avoids some computational problems near acos(1)
 theta = real ( acos(innerProd' ./ l)) ;  

end

