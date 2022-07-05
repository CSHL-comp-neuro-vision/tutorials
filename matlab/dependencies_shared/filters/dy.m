function result = dy(im,edges)
% DY: Y-derivative of an image, using upConv, and 5-tap derivative filters
%     designed by Eero Simoncelli.
% 
%      result=dy(im,edges)
% 
%      edges - any valid edge handler for upConv (default is 'circular').
%
% DJH '96

if ~exist('edges')
  edges='circular';
end

s = [3.342604e-2 0.241125 0.450898 0.241125 3.342604e-2];
sd = [9.186104e-2 0.307610 0.00000 -0.307610 -9.186104e-2];

tmp = upConv(im,sd',edges);
result = upConv(tmp,s,edges);
return;

%%% Debug
im=mkDisc(64);
res=dy(im);
displayImage(res)
