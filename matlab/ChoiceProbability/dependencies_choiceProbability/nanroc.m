function [a,p] = nanroc(x,y,N,min4roc)
% nanroc computes area under ROC from 2 distributions x and y
%    if x and y are matrices, nanroc returns a row vector of ROC vals for
%    each pair of column vectors.  nanroc removes any nan terms and makes
%    calls to rocN.  The N argument determines the number of points in the
%    ROC (default is 100). The min4roc is the number of finite values
%    necessary to compute the ROC (default is 3 per distribution). 
%	usage area = nanroc(x,y,N,min4roc)
%      
%     [area, p] = nanroc(x,y,N,min4roc,n4shuf)
%       will also return a pvalue for each roc area, testing against H0:
%       area=0.5.  The 5th argument n4shuf is the number of points in the
%       monte carlo distribution (permutation).  It is 100 by default.
%
%  

% update 1/10/96  added the plotflag option
% 1/24/98  added p vals from rocshuf calls (optionally)
% update 6/17/10 changed finite to isfinite

if nargin < 3
  N = 100;
  min4roc = 4;
end
if nargin < 5
  n4shuf = 100;				%  default for permutation test
end

% check matrix dimensions
[m n] = size(x);
[mm nn] = size(y);
if n~=nn
  error('nanroc: matrices must have the same number of columns')
end


a = nans(1,n);
p = nans(1,n);
% loop through columns
for i = 1:n
  x0 = x(:,i);
  L0 = isfinite(x0);
  x1 = y(:,i);
  L1 = isfinite(x1);
  if sum(L0) >= min4roc & sum(L1) >= min4roc
    if nargout < 2
      a(i) = rocN(x0(L0),x1(L1),N);
    else
      % compute area and pval
      [a(i) p(i)] = rocshuf(x0(L0),x1(L1),n4shuf);
    end
  end
end



