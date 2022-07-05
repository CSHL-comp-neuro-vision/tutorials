% this tutorial is a quick and dirty demonstration of how center-
% surround spatial receptive fields decorrelate images in space
% 
% things to try:
%  1. compare the results on different kinds of images (supplied)
%  2. try making a better receptive field model to get better results
%
% Pam Reinagel 2004

clear,pack,close all
current=pwd; % the present directory
path(path,pwd);

% get a sample image to play with
[x map]=imread('adams1.bmp','bmp');  
x=double(x); % required by conv2 later

% make a 15x15 square center surround filter
simpleRF=zeros(15); 
simpleRF(6:10,6:10)=ones(5);
simpleRF=simpleRF-mean(mean(simpleRF)); % set to mean 0, sum to 0
simpleRF=simpleRF*.5; % reduce amplitude

% version of image filtered by center/surround rf
procX=conv2(x-mean2(x),simpleRF); % it's padded with 7 pixels every side
procX=procX(8:size(x,1)+7, 8:size(x,2)+7);

figure
subplot(221), imshow(x,map)
subplot(222), imagesc(simpleRF), axis off
subplot(223), imagesc(procX), axis off


% make a large array for covaraince calcs
% this is a hack because of boundaries but good enough for illustration
% one for the raw images
x2=reshape(x',64,4800)'; % save all data as many rows of 64 pixels
x2=x2-mean2(x2); % remove DC

% and for the filtered (processed) versions
procx2=reshape(procX',64,4800)'; % save all data as many rows of 64 pixels
procx2=procx2-mean2(procx2); % remove DC

% compute correlations
Correlation=cov(x2); % covariance matrix for raw images
ProcCorrelation=cov(procx2); % covariance matrix for filtered images

% extract  square from images, required by OAP function
x=x(1:256,1:256); 
procX=procX(1:256,1:256);
procX=procX-min(min(procX)); % make procX positive
procX=procX/(max(max(procX))/255); %  normalize the range

%power spectrum
OAP= orientavg(x, 3);
ProcOAP=orientavg(procX, 3);
freqUnits = (1:128) * 0.2;% arbitrary units of degrees

figure
% show image
subplot(231); imshow(x,map); axis square, axis off
title('Natural Image')

% show correlations
subplot(232),
plot(-31.5:31.5,Correlation(32,:)/max(Correlation(32,:)))
axis square, box off
axis([-32 32 0 1]), hold on, plot([-32 32],[.5 .5],':');
xlabel('distance'), title('Correlation')

% the power spectrum
subplot(233),
loglog(freqUnits, OAP); axis square
title('Power spectrum')
xlabel('cycles/degree')

%now the processed image
subplot(234),hold off
imagesc(procX)
axis square, axis off
title('Filtered Image')

% show correlations
subplot(235),
plot(-31.5:31.5,ProcCorrelation(32,:)/max(ProcCorrelation(32,:)))
axis square, box off
axis([-32 32 0 1]), hold on, plot([-32 32],[.5 .5],':');

% the power spectrum
subplot(236),
loglog(freqUnits, ProcOAP); axis square
