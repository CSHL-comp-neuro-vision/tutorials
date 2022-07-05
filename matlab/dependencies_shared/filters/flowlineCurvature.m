function result = flowlineCurvature(im,edges)
% flowlineCurvature: One of two 2nd-order derivative-based
% measures that is invariant to any monotonic intensity
% nonlinearity.  From "Geometry-Driven Diffusion in Computer
% Vision", Romeny (ed.), p. 42.  See also isophoteCurvature.
%
%      result = flowlineCurvature(im,edges)
%      edges - any valid edge handler for upConv (default is 'circular').
%
% DJH 8/96

if ~exist('edges')
  edges='circular';
end

Ix=dx(im,edges);
Iy=dy(im,edges);
Ixx=dxx(im,edges);
Iyy=dyy(im,edges);
Ixy=dxy(im,edges);

result = ((Ix.^2 - Iy.^2).*Ixy - Ix.*Iy.*(Iyy - Ixx)) ./ ...
         ((Ix.^2 + Iy.^2).^(3/2)); 
     
return;

%%% Debug
ein=pgmRead('einstein.pgm');
res=flowlineCurvature(blur(ein,2));
displayImage(res,[-1,1])
