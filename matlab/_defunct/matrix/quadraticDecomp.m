function T = quadraticDecomp(M)
% QUADRATIC DECOMPOSITION: Given symmetric matrix M, want to get a matrix T
% so that T' M T = I.  To do this, write M = A A', take the SVD, and use
% what you get back.
%
% T = quadraticDecomp(M)
%
% DJH, 1/96.

[u,s,v] = svd(M);
T = u * sqrt(s);

