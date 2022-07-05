function x = invWeibull(w,a,b,type)
%
%	x = invWeibull(w,a,b,type)
%
%   Compute the inverse of a Weibull function
%     x:  the intensity levels
%     a,b:  the alpha and beta parameters in the weibull
%    type:  tafc or yesno. 
%           If tafc,  the form   1- 0.5*exp(- (x/a)^b)  is used.
%           if yesno, the form   1- exp(- (x/a)^b)      is used.
%     
% Modified from weibull.m by Xuemei Zhang

if nargin < 4
 type = 'tafc';
end

if type(1:4) == 'tafc'
  x = a * (-log((1-w)*2)).^(1/b);
elseif type(1:4) == 'yesn'
  x = a * (-log(1-w)).^(1/b);
end
