function result = sapiroGrad(f,edges,epsilon)
% SAPIROGRAD: Gradient function used by anisoSapiro.m 
% 1/(gradient+epsilon) (fy^2 fxx + fx^2 fyy - 2 fx fy fxy)^(1/3)
%
% result = sapiroGrad(im,edges,epsilon)
% im - input image
% edges- any vali edge handler for upConv (default is 'circular')
% epsilon - avoids divide by zero
%
% DJH '96

if ~exist('edges')
  edges='circular';
end

if ~exist('epsilon')
  epsilon = 1e-6;
end

fx = dx(f,edges);
fxx = dxx(f,edges);
fy = dy(f,edges);
fyy = dyy(f,edges);
fxy = dxy(f,edges);
fx2 = fx.^2;
fy2 = fy.^2;
grad = sqrt(fx2 + fy2);
g = 1 ./ (grad + epsilon);
tmp =fy2.*fxx + fx2.*fyy - 2.*fx.*fy.*fxy;
result = g.* sign(tmp) .* (abs(tmp)).^(1/3);

