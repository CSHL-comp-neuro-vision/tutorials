% anisoTest: test code for anisoPerona and anisoSapiro

% Noisy step
step = 1+[ones([32 64]) ; zeros([32 64])];
noise = 0.1 * randn(size(step));
im = (step + noise);
figure(1)
displayImage(im);
truesize

sap = anisoSapiro(im,40); 
figure(2)
displayImage(sap);
truesize

per = anisoPerona(im,20,0.1); 
figure(3)
displayImage(per);
truesize

figure(4)
clf
plot([im(:,32) sap(:,32) per(:,32)]);

figure(4)
clf
hist(sap(:))
hist(per(:))
