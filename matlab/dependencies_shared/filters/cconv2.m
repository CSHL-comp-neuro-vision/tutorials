function result = cconv2(im,filter)
% CCONV2: Circular convolution (calls upConv).
%
% res = cconv2(image,filter)
%
% DJH, 8/96
% update 12/97

result = upConv(im,filter,'circular',[1,1]);
return;

%%% Debug
im=mkImpulse(7);
filter = [1,2,4,2,1]'*[1,2,4,2,1];
filter=filter/sum(sum(filter));
res=cconv2(im,filter);
