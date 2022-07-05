function [C, Q] = qrancorrelmtx(n,rmin,rmax,nocheck)
% QRANCORRELMTX returns a random n by n correlation matrix, C, and its
%	matrix square root Q.  It does this by directly computing the elements
%	according to a formula derived for constant r value. 
%	This derivation is in my Mathematica file sqrt_matrix_normalization.math
%	The routine uses a uniform distribution of r values on rmin..rmax.
%	It goes through the actual matrix square root and is therefore slow unless
%	nocheck is set to nonzero.  Then the calculated matrices are returned.
%	Usage: [C Q] = rancorrelmtx(n,rmin,rmax,nocheck)

%	(M N Shadlen, 6/93)
% 7/23/97  mns forced symmetry on the random vals for correlation and hence
% on Q. 

rq = (rmax - rmin);
r = (rmin + rmax)/2;
won = ones(n);
R = (rmin * won) + (rq * rand(n));

% changes 7/23/97 mns

q = triu(R,1);
meanr = sum(q(:)) / (.5 * n*(n-1));
R = q + q' + meanr * eye(n);


% mean(R(:))
Q = sqrt(((2/n) * won) + R - (((2/n) * won) .* R) - ( 2 * sqrt(won - R) .* sqrt(won - R + n*won.*R)/n)) / sqrt(n);
Q = (~eye(n)) .* Q + eye(n) * sqrt(2/n + r - 2*r/n - 2*sqrt(1-r)*sqrt(1 - r + n*r)/n)* (2/r + 2*sqrt(1 - r)*sqrt(1 - r + n*r)/r)/(2*sqrt(n));
C = (~eye(n)).* (Q * Q') + eye(n);
% we can clean up by setting the diags to 1, taking square root and squaring again.
if nocheck ~= 1
  Q = sqrtm(C);
  C = Q * Q';
end


