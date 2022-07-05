function rm = meanR(x)
% rm = meanR(x)
%   given a matrix,x, calculate the correlation matrix and return the
%   average r value.  x is arranged so that the rows are observations and
%   the columns are the variables.  SO a 100 by 2 yields a simple 2 by 2
%   correlation matrix.
r = corrcoef(x);
n = length(r);
rm = (sum(r(:))-n)/(n*(n-1));

