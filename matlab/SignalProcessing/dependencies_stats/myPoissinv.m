function x = myPoissinv(p,lambda)
%myPoissinv Inverse of the Poisson cumulative distribution function (cdf).
%   X = POISSINV(P,LAMBDA) returns the inverse of the Poisson cdf 
%   with parameter lambda. Since the Poisson distribution is discrete,
%   POISSINV returns the smallest value of X, such that the poisson 
%   cdf evaluated, at X, equals or exceeds P.
%
%   The size of X is the common size of P and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
% Modified from Matlab poissinv (DJH, 9/97) to fix a bug (see
% poissinv.m for example of the bug).

% Bug in Mathworks poissinv
%   poissinv(0.3,1)
%   poissinv(0.4,1)
%   poissinv(0.8,1)
%   poissinv([0.3,0.4,0.8],1)
%
%   myPoissinv(0.3,1)
%   myPoissinv(0.4,1)
%   myPoissinv(0.8,1)
%   myPoissinv([0.3,0.4,0.8],1)


if nargin < 2, 
    error('Requires two input arguments.'); 
end

[errorcode p lambda] = distchck(2,p,lambda);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

x = zeros(size(p));

count = 0;
cumdist = poisspdf(count,lambda);

% Compare P to the poisson cdf.
k = find(lambda > 0 & p >= 0 & p < 1);
while any(any(p(k) > cumdist(k)))
    index = find(cumdist(k) < p(k));
    x(k(index)) = x(k(index)) + 1;  
    count = count + 1;
    cumdist(k) = cumdist(k) + poisspdf(count,lambda(k));
end

% Return NaN if the arguments are outside their respective limits.
k = find(lambda <= 0 | p < 0 | p > 1);
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k)));
end

% Return Inf if p = 1 and lambda is positive.
k = find(lambda > 0 & p == 1);
if any(k)
    tmp  = Inf;
    x(k) = tmp(ones(size(k)));
end

