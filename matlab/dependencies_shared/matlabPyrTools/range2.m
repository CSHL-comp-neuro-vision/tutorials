function [MIN,MAX] = range2(MTX)
% [MIN, MAX] = range2(MTX)
%
% Compute minimum and maximum values of MTX, returning them as a 2-vector.

% Eero Simoncelli, 3/97 as a mex function
% Converted to matlab function 7/2019
MIN = min(MTX(:));
MAX = max(MTX(:));
