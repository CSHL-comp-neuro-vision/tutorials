function[ warp_im ] = warpPolar( im, theta, mag, method )
% WARPPOLAR
%
% USAGE: warpPolar( im, theta, mag, method );
%
% DESCRIPTION: The function warp takes an image, IM, an angle and magnitude
%		image, THETA and MAG respectively, of the same size as IM,
%		and an interpolation method, METHOD = {'nearest', 'bilinear',
%		'cubic'}.  The function warps IM according to the specified
%		displacement field (i.e., THETA and MAG).
%
% NOTES: The built-in interp2() routine places NaN in the warped image
%	  where it is unable to interpolate.  For display purposes, this
%	  routine replaces the NaN's with 0's.
%
% EPS '96


  if (exist('method') ~= 1)
		method = 'bicubic';
	end
		
	[ydim, xdim] 	= size(im);
	[x_ramp,y_ramp]	= meshgrid( 1:xdim, 1:ydim );
	warp_x		= x_ramp + (mag .* cos( theta ));
	warp_y		= y_ramp - (mag .* sin( theta ));
	warp_im		= interp2(warp_x, warp_y, im, x_ramp, y_ramp, method);

	warp_im		= reshape( warp_im, 1, ydim*xdim );
	nan_i		= find( isnan( warp_im ) );
	warp_im( nan_i )= zeros(1, size(nan_i, 2) );
	warp_im		= reshape( warp_im, ydim, xdim );

	display_im( warp_im );

