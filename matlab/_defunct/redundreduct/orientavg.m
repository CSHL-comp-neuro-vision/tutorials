% function [OrientationAverage] = f(imblock, normalize, plotresult, sampling);
%
%  inputs:
%    imblock    is a square image block with a linear intensity map
%    normalize  if 0 don't normalize 
%               if 1 only subtract the mean intensity (remove DC) DEFAULT
%               if 2 also divide by mean intensity (normalize to contrast=total power)
%               if 3 do exactly as Hans vH. (subtract *weighted* mean not DC)
%
%    plotresult is 0 if no plot, 1 if plot, default 0
%    sampling   is the spatial sampling frequency in cycles/degree
%               sampling is only used to label plot axis
%               defaults to global SamplingDefault
%  output:
%    OrientationAverage is a vector containing the average over
%         orientations of the PSD at each of BlockSize/2 freqs.
%         the PSD is calculated using a conical window

function [OrientationAverage] = f(imblock, normalize, plotresult, sampling);

global SamplingDefault;

if(nargin<2) normalize=1; end
if(nargin<3) plotresult=0; end
if(nargin<4) sampling=SamplingDefault; end

BlockSize = size(imblock,1);

%%%%%%%%% set up %%%%%%%%%%%%%%

%map each x,y pair to a linear list of pixels
pixels = (1:BlockSize^2);
x = ceil(pixels/BlockSize);
y = pixels-(x-1)*BlockSize;

%calculate dist from center of grid to center of pixel
%  BlockSize/2,BlockSize/2 is the center of the grid,
%  x-.5, y-.5 is the center of the pixel with coords x,y
dist=sqrt( (BlockSize/2 - (x - .5) ).^2 + ...
         (BlockSize/2 - (y - .5) ).^2 );

%%%%%% make window %%%%%
window=zeros(BlockSize);
for i=1:BlockSize^2, 
  window(x(i),y(i)) = dist(i);                       % make square array from linear
end 
window=window/ (BlockSize/2);                        % define dist to edge as 1
window=ones(BlockSize) - window;                     % now center is 1, edge is 0
window(find(window<0))=zeros(size(find(window<0)));  % set corners to 0 not negative


if (normalize>0) % normalize for DC
  if (normalize < 3) % just subtract mean
    mu_imblock = mean2(imblock);
  else % use Hans's weighted mean
    mu_imblock = sum(sum( imblock .* window)) / sum(sum(window));
  end
  imblock = imblock - ones(size(imblock)) * mu_imblock;
end


if (normalize>1) %normalize for total contrast
  imblock = imblock ./ mu_imblock;
end

%%%% calc the psd %%%%
PSD = abs(fft2(imblock .* window)).^2 / (norm(window))^2;
Centered  = fftshift(PSD);

%%%% orientation average %%%%
for i=1:floor(BlockSize/2), % as many freq bins as in the fft
     clear indx indy;         % because not same size each time
     clear FreqBandPts;

     indx = x(find(dist>i-1 & dist<=i));
     indy = y(find(dist>i-1 & dist<=i));
     
    FreqBandPts = diag(Centered(indx,indy));
    OrientationAverage(i) = sum(FreqBandPts)/size(FreqBandPts,1);

end

if (plotresult > 0),
   freqUnits = (1:floor(BlockSize/2)) * sampling;
   h=figure; plot(log10(freqUnits),log10(OrientationAverage)), 
   title('Orientation Average of PSD'), ylabel('log10(Power)')
   xlabel('Log Spatial Frequency (cycles/degree)'), figure(h);
end



