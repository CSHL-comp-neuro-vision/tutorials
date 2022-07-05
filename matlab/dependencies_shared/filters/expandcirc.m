% EXPANDCIRC	image expansion followed by circular convolution
%
%	expandcirc(im,filter) performs a 1D or 2D circular convolution
%	of image im with filter.
%
%	expandcirc(im,filter,[stepx,stepy],[startx,starty],[s_x,s_y]) performs
%	convolution with expansion step sizes (stepx,stepy), starting
%	point (startx,starty) in the image, and output image size
%	(s_x,s_y).  
%
% Provided for backward compatibility.  Calls upConv.
%
% DJH, 12/97

function result = expandcirc(im,filter,step,start)

if (exist('step') ~= 1)
  step = [1,1];
end	

if (exist('start') ~= 1)
  start = [1,1];
end	

result = upConv(im,filter,'circular',step,start);
return;

%%% Debug
im=mkImpulse(7);
filter = [1,2,4,2,1]'*[1,2,4,2,1];
filter=filter/sum(sum(filter));
res=expandcirc(im,filter);
res=expandcirc(im,filter,[2 2]);
