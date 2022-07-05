function result = anisoSapiro(im,numiters,edges,deltaT,epsilon)
% anisoSapiro: Anisotropic diffusion of an image, using Guillermo Sapiro's
%     method.
%
% result = anisoSapiro(im,numiters,edges,deltaT,epsilon)
% 
%      result=anisoSapiro(im,numiters,deltaT,epsilon)
%      im - input image
%      numiters - number of iterations
%      edges - any valid edge handler for upConv (default is 'circular').
%      deltaT - time step (default 0.02)
%      epsilon - avoids divide by zero in sapiroGrad (default 1e-6)
%
% Example in anisoTest.m
%
% DJH '96

if ~exist('deltaT')
  deltaT = 0.02;
end

if ~exist('edges')
  edges='circular';
end

if ~exist('epsilon')
  epsilon = 1e-6;
end

result = im(:,:);
for iter = 1:numiters
  result = result + deltaT .* sapiroGrad(result,edges,epsilon);
end
