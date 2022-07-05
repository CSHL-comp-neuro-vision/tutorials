% test.m
%
% Tests the function CONVOLVECIRC, UPCONV and CORRDN by convolving a random
% noise image with a low-pass filter.
%
% gmb, 6/1/96
% djh, 8/1/96, modified to include upconv and corrdn

%create a random noise image
im=round(rand(100,100));

%show it
figure(1)
imagesc(im);
colormap(gray);
truesize;
axis off

%create a filter.
filter = [1,2,4,2,1]'*[1,2,4,2,1];
filter=filter/sum(sum(filter));

%filter the image
foo=convolvecirc(im,filter);
bar=expandcirc(im,filter);

%show the filtered image
imagesc(foo);
imagesc(bar);

%%-----------------------------------------------------------------

im = eye(7);
im = eye(8);

corrDn(im,2.^[0:3]',[1,1])
upConv(im,2.^[0:3]',[1,1])

%%timing:
im = rand(100); filt = rand(20);
tic; corrDn(im,filt); toc
tic; corrDn(im,filt,[1,1],[0,0],'reflect1'); toc
tic; rconv2(im,filt); toc
tic; conv2(im,filt); toc
