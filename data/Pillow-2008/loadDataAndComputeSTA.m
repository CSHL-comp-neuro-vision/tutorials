% load and compute spike-triggered average of one example neuron

load Stim_reduced;
load SpTimesRGC;

cellnum = 1; % which cell to examine
nlags = 20; % number of time bins of lags to consider

nx = 10; % number of spatial pixels in x direction
ny = 10; % number of spatial pixels in y direction
npx = nx*ny;

slen = size(Stim,1); % number of time bins

% bin spike times
sps = hist(SpTimes{1},0.5:1:slen)';

% compute STA
sta = zeros(nlags,npx);
for jj = 0:nlags-1
    sta(nlags-jj,:) = sps(jj+1:end)'*Stim(1:end-jj,:);
end
sta = sta/sum(sps); % divide by # spikes

subplot(221);
imagesc(sta);
xlabel('space'); ylabel('time');
title('full STA');

subplot(223)
plot(sta(:,36));
xlabel('time'); 
title('time slice (pixel 36)'); 

subplot(224)
imagesc(reshape(sta(17,:),nx,ny));
xlabel('x'); ylabel('y');
title('space slice (lag 3)');

