function [outimage] = anisoPerona(inimage,iterations,K)
%
% anisoPerona: From the Chapter "Anisotropic Diffusion", by Perona, Shiota,
% and Malik, in Geometry-Driven Diffusuion, Bart M. ter Haar Romeny, p. 81.
%
% PURPOSE: Perform a smoothing over a gray scale image
%
% [outimage] = anisoPerona(inimage,iterations,K)
% 
% Example in anisoTest.m


if nargin < 3
 K = 2;
end
if nargin < 2
 iterations = 8;
end

lambda = 0.25;
outimage = inimage;       [m,n] = size(inimage);

rowC = [1:m];         rowN = [1 1:m-1];    rowS = [2:m m];
colC = [1:n];         colE = [1 1:n-1];    colW = [2:n n];

for i = 1:iterations,

  deltaN = outimage(rowN,colC) - outimage(rowC,colC);
  deltaE = outimage(rowC,colE) - outimage(rowC,colC);

  fluxN = deltaN .* exp( - (1/K) * abs(deltaN) );
  fluxE = deltaE .* exp( - (1/K) * abs(deltaE) );

  outimage = outimage + lambda * ...
    (fluxN - fluxN(rowS,colC) + fluxE - fluxE(rowC,colW));

end;


