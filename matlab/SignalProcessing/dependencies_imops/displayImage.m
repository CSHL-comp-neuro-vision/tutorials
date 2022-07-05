function retval = displayImage( im, range, nshades, style );
% DISPLAYIMAGE - Display matlab matrices as images.
%
% range = displayImage (mtx, range, nshades, style )
% range = displayImage (R,G,B, style )
% range = displayImage (Y,I,Q, style )
% 
% Display a MatLab matrix in the current figure, as a grayscale image.
% 
% If the matrix is complex, the real and imaginary parts are shown
% side-by-side.
% 
% RANGE (optional) is a 2-vector specifying the values that map to
% black and white, respectively. Passing a value of 'AUTO1' (default)
% auto-scales to fill range.  'AUTO2' puts image mean at 50% gray.
% 
% NSHADES (optional) specifies the number of gray shades, and defaults
% to the size of the current colormap.
%
% STYLE (optional) defaults to 'gray' which is the behavior above.
% Other settings are 'rgb' and 'yiq' which allow color images to
% be displayed. The color channels are given as seperate matrix parameters.
% If three fullsize matrices are given, 'rgb' will be assumed. 
%
% Also multiple images can be displayed in a figure with subplot(), and
% using the (m)ulti styles. 'mgray', 'mgrb', 'myiq'. These styles simply
% keep displayImage from reseting the figure for you.
% 
% See also RGB2NTSC, SUBPLOT, TIFFREAD
%
% JED 12/96 : EPS 6/96 


if (exist('range') ~= 1)
  range = 'auto1';
end

if (exist('nshades') ~= 1)
  nshades = size(colormap,1);
  if (nshades>1024)
	nshades=256;
  end;
end

% Handle optional extra display modes
if (exist('style') ~= 1)
	if(prod(size(im))==prod(size(nshades)))
		style='rgb';
	else
		style='gray';
	end;	
end;

if (style(1:3)=='mgr' | style(1:3)=='gra')
	if (nshades < 2)
	  nshades = 2;
	end
end


if (style(1:3)=='rgb' | style(1:3)=='mrg')
	[d_im map]=rgb2ind(im,range,nshades);
	colormap(map);
	retval(1)=min(min(im));
	retval(2)=max(max(im));
	retval(3)=min(min(range));
	retval(4)=max(max(range));
	retval(5)=min(min(nshades));
	retval(6)=max(max(nshades));
end;

if (style(1:3)=='yiq' | style(1:3)=='myi')
	sz = prod(size(im));
	YIQ = [reshape(im,sz,1),reshape(range,sz,1),reshape(nshades,sz,1)];
	RGB = ntsc2rgb(YIQ); 
	clear YIQ;
	R = reshape(RGB(:,1),size(im));
	G = reshape(RGB(:,2),size(im));
	B = reshape(RGB(:,3),size(im));
	clear RGB;
	[d_im map]=rgb2ind(R,G,B);
	clear R; clear G; clear B;
	colormap(map);
end;

if (style(1:3)=='gra' | style(1:3)=='mgr')
	if (strcmp(range,'auto1') | strcmp(range,'auto'))
	  if isreal(im)
	    mn = min(min(im));
	    mx = max(max(im));
	  else
	    mn = min(min(min(real(im))),min(min(imag(im))));
	    mx = max(max(max(real(im))),max(max(imag(im))));
	  end
	  range = [mn,mx];
	elseif strcmp(range,'auto2')
	  if isreal(im)
	    mn = min(min(im));
	    mx = max(max(im));
	    av = mean2(im);
	  else
	    mn = min(min(min(real(im),imag(im))));
	    mx = max(max(max(real(im),imag(im))));
	    av = (mean2(real(im)) + mean2(imag(im)))/2;
	  end
	  mx = max(mx-av,av-mn);
	  range = [av-mx,av+mx];
	end    

	if ((range(2) - range(1)) <= eps)
	  range(1) = range(1) - 0.5;
	  range(2) = range(2) + 0.5;
	end

	colormap(gray(nshades));

	if isreal(im)
	  factor=1;
	else
	  factor = 1+sqrt(-1);
	end

	%% same as (im-mn)/(mx-mn) + 1:
	d_im = im * ((nshades-1) / (range(2)-range(1))) + ...
	    factor*(1.5 - (range(1)*(nshades-1))/(range(2)-range(1)));
	
	retval=range;
end;


if isreal(d_im)
  if (style(1:1)~='m'), subplot(1,1,1); end;
  hh = image( d_im );
  if (style(1:3)=='gra' | style(1:3)=='mgr')
    set(hh,'UserData',range);
  end;
  axis('image');
  axis('off');
else
  subplot(1,2,1);
  hh = image(real(d_im));
  set(hh,'UserData',range);
  axis('image');
  axis('off');
  subplot(1,2,2);
  hh = image(imag(d_im));
  set(hh,'UserData',range);
  axis('image');
  axis('off');
end  


