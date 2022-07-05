function [out,xFilt,thetaFilt] = convolveImage(p,img,xWidth,thetaWidth)
%[out,xFilt,thetaFilt] = convolveImage(p,img,xWidth,thetaWidth)

dx = p.x(2)-p.x(1);
dtheta = p.theta(2)-p.theta(1);

%standard convolution for space (x) dimension:
xc = p.x(ceil(length(p.x)/2));
xFilt =  normpdf(p.x,xc,xWidth);
out = dx*conv2(img,xFilt,'same');

%circular convolution for orientation (theta) dimension
thetac = p.theta(ceil(length(p.theta)/2));
thetaFilt = normpdf(p.theta,thetac,thetaWidth);

thetaFilt = repmat(thetaFilt(:),1,length(p.x));
out = dtheta*fftshift(real( ifft( fft(out).*fft(thetaFilt))),1);