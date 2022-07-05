function [mn,mx] = range2(X)
%[MIN, MAX] = range2(X)
% Compute minimum and maximum values of matrix X, returning them as a 2-vector.

mn = min(X(:));
mx = max(X(:));

