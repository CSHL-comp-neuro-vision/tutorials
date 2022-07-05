function [newImage, newSize] = magnifyImage(im, imSize, magFactor )
%
%   [newImage, newSize] = magnifyImage(im, imSize, magFactor )
%
%PURPOSE:
%   Magnify an image to a larger size smoothly using int
%   bilinear interpolation (interp2).  The routine assumes that
%   the image is encoded using a standard gray map.  General indexed
%   images must be converted to this format for this to make sense.
%
%      im:  an image in the mr format, (long vector)
%      imSize: the image x and y dimensions
%      magFactor:  The approximate magnification factor for the image
%
%AUTHOR:
%   Denise Krol
%
%DATE:  08.15.94
%

x = [1:imSize(2)];
y = [1:imSize(1)]';

newSize = magFactor .* imSize;

xi = linspace(1,imSize(2), newSize(2));
yi = linspace(1,imSize(1), newSize(1))';
z = reshape(im,imSize(1),imSize(2));

newImage = interp2(x,y,z,xi,yi);
newSize = size(newImage);
newImage = reshape(newImage,1,newSize(1)*newSize(2));
























