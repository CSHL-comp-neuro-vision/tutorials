function result = dyy(im,edges)
% DXX: 2nd Y-derivative of an image, using upConv, and 5-tap derivative
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

s = [3.342604e-2 0.241125 0.450898 0.241125 3.342604e-2];
sdd = [0.202183 9.181186e-2 -0.587989 9.181186e-2 0.202183];

tmp = upConv(im,sdd',edges);
result = upConv(tmp,s,edges);
return;

%%% Debug
im=mkDisc(64);
res=dyy(im);
displayImage(res)
