function x = expinv(p,mu)
%EXPINV Inverse of the exponential cumulative distribution function.
%   X = EXPINV(P,MU) returns the inverse of the exponential 
%   cumulative distribution function, with parameter MU,
%   at the values in P.
%
%   The size of X is the common size of P and MU. A scalar input   
%   functions as a constant matrix of the same size as the other input.    

%     Copyright 1993-2000 The MathWorks, Inc. 
%     $Revision: 2.9 $  $Date: 2000/05/26 17:28:52 $


% Set the default mean to 1.
if nargin < 2,
    mu = 1;
end

if nargin < 1, 
    error('Requires at least one input argument.'); 
end

[errorcode p mu] = distchck(2,p,mu);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
x=zeros(size(p));

%   Return NaN if the arguments are outside their respective limits.
k = find(p < 0 | p > 1 | mu <= 0);
if any(k),
   tmp  = NaN; 
   x(k) = tmp(ones(size(k)));
end

k = find(p > 0 & p < 1);
if any(k),
    x(k) = -mu(k) .* log(1 - p(k));
end
