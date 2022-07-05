% RES = warpIm( IM, Ywarp, Xwarp, METHOD );
%
% Warp image IM according to vector field contained in [Ywarp,Xwarp]:
%    RES(y,x) = IM(y+ywarp,x+xwarp)
%
% METHOD determines the interpolation method.  Choices are
%     {'nearest', 'linear', 'cubic'}, default='cubic'
%
% NOTES: The built-in interp2() routine places NaN in the warped image
%	  where it is unable to interpolate.  For display purposes, this
%	  routine replaces the NaN's with 0's.

% Hany Farid, Spring 96.
% Modified by EPS, 8/97: Arg order is now (y,x), and does a "from" warp.

function [ res ] = warp( im, ywarp, xwarp, method )

if (exist('method') ~= 1)
  method = 'cubic';
end
		
[ydim, xdim] 	= size(im);
[x_ramp,y_ramp]	= meshgrid( 1:xdim, 1:ydim );
warp_x		= x_ramp - xwarp;
warp_y		= y_ramp - ywarp;
res		= interp2(x_ramp, y_ramp, im, warp_x, warp_y, method);

res		= reshape( res, 1, ydim*xdim );
nan_i		= find( isnan( res ) );
res( nan_i )= zeros(1, size(nan_i, 2) );
res		= reshape( res, ydim, xdim );

