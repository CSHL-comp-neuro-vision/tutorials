function z = multiRandn(MeanVec,Cov)
% MULTIRANDN: Draws a vector from a multivariate normal distribution, given
% the mean and covariance.
%
% z = multiRandn(MeanVec,Cov)
%
% DJH, 1/96.

MeanVecCol = MeanVec(:);
A = quadraticDecomp(Cov);
y = randn(size(MeanVecCol));
y = A * y + MeanVecCol;
z = reshape(y,size(MeanVec,1),size(MeanVec,2));
