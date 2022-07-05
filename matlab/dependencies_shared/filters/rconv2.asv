function result = rconv2(im,filt)
% RCONV2: Convolution, with boundaries handled via
% reflection about the edge pixels (calls upConv).
%
% result = rconv2(image,filter)
%
% DJH, 8/96

result = upConv(im,filt,'reflect1');


%%% Debug
% im=mkImpulse(7);
% filter = [1,2,4,2,1]'*[1,2,4,2,1];
% filter=filter/sum(sum(filter));
% res=rconv2(im,filter);

