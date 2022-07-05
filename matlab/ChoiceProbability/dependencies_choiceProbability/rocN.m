function a = rocN(x,y,N,plotflag)
% ROC	computes area under ROC given distributions x and y
%	uses N points to construct the ROC. Default is 100. 
%	if plotflag ==1 the roc curve is drawn in the current figure. The
%	default is *not* to plot the curve.
%	usage area = roc(x,y,N,plotflag)
%    Note: x and y must contain at least 3 values.  Function returns nan if
%    insufficient data are supplied.

% update 1/10/96  added the plotflag option
% 1/21/97  return nan for empty inputs and vector length < 3

if nargin < 3
  N = 100
end
if isempty(x) | isempty(y)
  a = nan;
  return;
end
[m n] = size(x);
x = reshape(x,1,m*n);
[m n] = size(y);
y = reshape(y,1,m*n);

if length(x)<3 | length(y)<3
  a = nan;
  fprintf(2,'Warning (rocN): too few trials. Function returns NaN');
  return;
end
zlo = min([min(x(:)) min(y(:))]);
zhi = max([max(x(:)) max(y(:))]);
z = linspace(zlo,zhi,N);
fa = zeros(1,N);	% allocate the vector
hit = zeros(1,N);
for i = 1:N
  fa(N-i+1) = sum(y > z(i));
  hit(N-i+1) = sum(x > z(i));
end
[m,ny] = size(y);
fa = fa/ny;
[m,nx] = size(x);
hit = hit/nx;
fa(1) = 0;
hit(1) = 0;
fa(N) = 1;
hit(N) = 1;
a = trapz(fa,hit);
% uncomment next line if you want to see the plot
if nargin > 3
  if plotflag == 1
    plot(fa,hit),axis('square'),xlabel('FA'),ylabel('Hit');
    title(sprintf('ROC area = %.3f', a));
  end
end



