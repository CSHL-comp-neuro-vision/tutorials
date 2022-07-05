function w = weibull(x,a,b,type)
%
%	w = weibull(x,a,b,type)
%
%AUTHOR:  Wandell
%DATE:  5.5.95
%PURPOSE:
%   Compute a Weibull function
%     x:  the intensity levels
%     a,b:  the alpha and beta parameters in the weibull
%    type:  tafc or yesno.  
%		          If tafc    the form   1- 0.5* exp(- (x/a)^b) is used.  (DEFAULT)help weib
%              if yesno, the form   1- exp(- (x/a)^b)         is used.
%     

if nargin < 4
 type = 'tafc';
end

if type == 'tafc'
  w = 1 - 0.5*exp( -(x/a).^b);
elseif type == 'yesno'
  w = 1 - exp( -(x/a).^b);
end
