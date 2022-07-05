function result = blur(im,levels,Filt)
% BLUR: Blurs an image by blurring and subsampling repeatedly, followed by
% upsampling and blurring.
%
%      result=blur(im,[levels],[Filt])
%
%      im - input image.
%      levels - number of times to blur and subsample (default is 1).
%      Filt - blurring 1d filter to be applied separably to the rows and
%               cols of im (default ='binom5').
%
% DJH '96
% update 12/97 to conform to Eero's updated pyrTools

if nargin < 2, levels = 1; end
if nargin < 3, Filt = 'binom5'; end
% if ~exist('levels')
%   levels=1;
% end
% 
% if ~exist('Filt')
%   Filt = 'binom5';
% end

if isstr(Filt)
  Filt = namedFilter(Filt);
end  

tmp = blurDn(im,levels,Filt);
result = upBlur(tmp,levels,Filt);

