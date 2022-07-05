function [ent, rmse, res] = qmfCompress(im, binsize, levels)
% qmfCompress: QMF image compression
%
% EPS '96
  
ysz = size(im,1);
xsz = size(im,2);
sz = xsz*ysz;
filt_sz = 9;

%% If levels argument was not included, compute it automatically:
if (~exist('levels'))
  levels = floor(log2(min(ysz,xsz)) - log2(filt_sz))
end

if (levels <= 0)  % low frequency band:

  res = im;  % don't compress
  ent = entropy2(res);

else

  lo = [ 0.02807382 -0.060944743 -0.073386624 0.41472545 0.7973934  ...
    0.41472545 -0.073386624 -0.060944743 0.02807382 ];
  hi = [ 0.02807382 0.060944743 -0.073386624 -0.41472545 0.7973934  ...
    -0.41472545 -0.073386624 0.060944743 0.02807382 ];

%% compute vertical, horizontal and diagaonal images.  Quantize using binsize.
  v_im = filter2(hi',filter2(lo,im));
  v_im = v_im(2:2:ysz,1:2:xsz);
  v_im = binsize * round(v_im/binsize);
  
  h_im = filter2(lo',filter2(hi,im));
  h_im = h_im(1:2:ysz,2:2:xsz);
  h_im = binsize * round(h_im/binsize);

  d_im = filter2(hi',filter2(hi,im));
  d_im = d_im(2:2:ysz,2:2:xsz);
  d_im = binsize * round(d_im/binsize);

%% compute lowpass image.  
  lo_im = filter2(lo',filter2(lo,im));
  lo_im = lo_im(1:2:ysz,1:2:xsz);

%% Call qmf_compress recursively on the lo_im.
  [lo_ent, lo_mse, lo_im] = qmfCompress(lo_im, binsize/2, levels-1);
  
%% reconstruct im from quantized low, vert, hor, diag images
  pred_im = zeros(ysz,xsz);
  pred_im(1:2:ysz,1:2:xsz) = lo_im;
  pred_im = filter2(lo',filter2(lo,pred_im));
  res = pred_im;

  pred_im = zeros(ysz,xsz);
  pred_im(2:2:ysz,1:2:xsz) = v_im;
  pred_im = filter2(hi',filter2(lo,pred_im));
  res = res + pred_im;

  pred_im = zeros(ysz,xsz);
  pred_im(1:2:ysz,2:2:xsz) = h_im;
  pred_im = filter2(lo',filter2(hi,pred_im));
  res = res + pred_im;

  pred_im = zeros(ysz,xsz);
  pred_im(2:2:ysz,2:2:xsz) = d_im;
  pred_im = filter2(hi',filter2(hi,pred_im));
  res = res + pred_im;

%% compute entropy: average of low, vert, hor, diag entropies
  ent = (lo_ent*length(lo_im(:)) + entropy2(v_im)*length(v_im(:)) + ...
      entropy2(h_im)*length(h_im(:)) + entropy2(d_im)*length(d_im(:)))/sz;
end

%% compute rmse:
rmse = sqrt(mean( (res(:) - im(:)).^2 ));
return;

%%% DISPLAY results:
if (levels > 0)
  max_levs = 5;
  ent/log2(var(v_im(:)))
  subplot(3,max_levs,levels); 
  imagesc(v_im); axis('image');
  subplot(3,max_levs,max_levs+levels);
  imagesc(h_im); axis('image');
  subplot(3,max_levs,2*max_levs+levels);
  imagesc(d_im); axis('image');
end

%%% Debug:
ein=pgmRead('einstein.pgm');
[ent,rmse,res]=qmfCompress(ein,10);
displayImage(ein+j*res,[0,255])
