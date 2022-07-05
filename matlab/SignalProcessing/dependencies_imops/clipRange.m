function [res, nClip]= clipRange(in,mn,mx)
%
%     res = clipRange(in,mn,mx)
%
%AUTHOR:  Wandell
%DATE:  08.16.95
%PURPOSE:
% Clip values in the input data to a minimum and maximum range.
% if mn and mx are not passed they are set to 0 and 1, respectively.
% 
%ARGUMENTS
%
%RETURNS
%

res = in;

if nargin < 3
 mx = 1;
end
if nargin < 2
 mn = 0;
end

l = find(in < mn);
nClip(1) = length(l);
res(l) = ones(1,nClip(1))*mn;

l = find(in > mx);
nClip(2) = length(l);
res(l) = ones(1,nClip(2))*mx;
