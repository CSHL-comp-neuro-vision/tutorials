function result = dxy(im,edges)
% DXY: Cross-derivative of an image, using upConv, and 5-tap derivative
%      filters designed by Eero Simoncelli.
% 
%      result=dx(im,edges)
% 
%      edges - any valid edge handler for upConv (default is 'circular').
%
% DJH '96

if ~exist('edges')
  edges='circular';
end

sd = [-9.186104e-2 -0.307610 0.00000 0.307610 9.186104e-2];

tmp = upConv(im,sd,edges);
result = upConv(tmp,sd',edges);
return;

%%% Debug
im=mkDisc(64);
res=dxy(im);
displayImage(res)
