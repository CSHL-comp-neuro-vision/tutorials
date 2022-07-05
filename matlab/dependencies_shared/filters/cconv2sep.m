function result = cconv2sep(im,rowfilt,colfilt)
% CCONV2SEP: Separable, circular, convolution using convolvecirc.
% 
%      result=cconv2sep(im,rowfilt,colfilt)
%
%      im - input image.
%      rowfilt - 1d filter applied to the rows
%      colfilt - 1d filter applied to the cols
%
% Example: foo=cconv2sep(im,[1 4 6 4 1],[-1 0 1]);
%
% DJH '96

rowfilt=rowfilt(:)';
colfilt=colfilt(:);

tmp = upConv(im,rowfilt,'circular',[1,1]);
result = upConv(tmp,colfilt,'circular',[1,1]);
return;

%%% Debug
im=mkImpulse(7);
filter = [1,2,4,2,1];
filter=filter/sum(filter);

res1=cconv2sep(im,filter,filter);
res2=cconv2(im,filter'*filter);
mse(res1,res2)
