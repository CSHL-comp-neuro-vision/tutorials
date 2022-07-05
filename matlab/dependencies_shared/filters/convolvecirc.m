% CONVOLVECIRC	circular convolution
%
%	convolvecirc(im,filter) performs a 1D or 2D circular convolution
%	of image im with filter.
%
%	convolvecirc(im,filter,[stepx,stepy],[startx,starty]) performs
%	convolution with subsampling step sizes (stepx,stepy) and starting
%	point (startx,starty) in the image.
%
% Provided for backward compatibility.  Calls corrDn.
%
% DJH, 12/97

function result = convolvecirc(im,filter,step,start)

if (exist('step') ~= 1)
  step = [1,1];
end	

if (exist('start') ~= 1)
  start = [1,1];
end	

result = corrDn(im,filter,'circular',step,start);
return;

%%% Debug
im=mkImpulse(10);
filter = [1,2,4,2,1]'*[1,2,4,2,1];
filter=filter/sum(sum(filter));
res=convolvecirc(im,filter);
res=convolvecirc(im,filter,[2 2]);
