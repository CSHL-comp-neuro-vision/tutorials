% -------------------------------------------------------------------------
% STCovTutorial.m

addpath('dependencies_STCovTutorial');
addpath('dependencies_motion');

% This tutorial demonstrates how spike-triggered covariance can
% be used to recover functional models in white-noise experiments.
% The technique is applied to the Adelson-Bergen opponent Energy 
% model. Before running this tutorial, you should run the 
% motionTutorial to aquaint yourself with spatiotemporal 
% receptive fields.
%
% For a more complete description of the method, see: Simoncelli et al
% (2004) Characterization of neural responses with stochastic stimuli.
% In: The New Cognitive Sciences, 3rd edition (ed. Gazzaniga). MIT Press.
%
% Nicole Rust, June 2004
% -------------------------------------------------------------------------


%% Gaussian white noise stimulus:
% The white noise stimlus we will be using to recover models of 
% simulated neurons is a bar stimulus.  The intensity of the bar
% values are chosen from a Gaussian distribution.
nT=11;
nX = 15;
stim = randn(100000,nX);
stim = stim./max(max(abs(stim)));

% Plot the stimulus presented on each frame and the stimulus history over
% the last 11 frames.  Notice that the random stimulus produces motion
% that is sometimes net leftward, sometimes rightward.
figure
for n=nT:nT+10
    subplot(2,2,1)
    %n
    imagesc(stim(n,:),[-1 1]); 
    colormap(gray);
    title('Stimulus on this frame')
    xlabel('space (bar number)')
    set(gca,'YTickLabel',[]);
    subplot(1,2,2)
    imagesc(flipud(stim(n-nT+1:n,:)),[-1 1]);
    title(['Stimuli presented on the last ',num2str(nT),' frames'])
    set(gca,'YTickLabel','0|-1|-2|-3|-4|-5|-6|-7|-8');
    xlabel('space (bar number)')
    ylabel('time (frames)')
    pause(0.25)
end;

% Our model cell will depend on the bar history on the 11 frames preceeding
% a spike. To model the cell, we begin by constructing a matrix whose rows 
% contain the 11-frame bar history at successive points in time:
allStim = zeros(size(stim,1)-nT+1, nX*nT);
for n = 1:size(allStim,1)
  allStim(n, :) = reshape(stim(n:n+nT-1,:),[1,nX*nT]);
end

% OPPONENT ENERGY MODEL --------------------------------------------------
%
% The energy model produces a phase-invariant response by combining the
% summed squared output of two filters in quaderature.  The opponent energy
% model computes the difference between a pair tuned for rightward motion
% and a second tuned for leftward motion.  Spike-triggered covariance
% recovers a model of the cell containing four equivalent filters.

% The linear filters used for the models in this tutorial are
% space-time oriented receptive fields.  For a description
% of space-time representations of motion and construction
% of directionally tuned linear filters, see the motionTutorial.
% We will be working with four filters: a pair tuned for rightward
% motion and a second pair tuned for leftward motion.  The members
% of each pair have a quaderature (90 degree) phase relationship:
[leftward_1, leftward_2, rightward_1, rightward_2] = ABmodel;

% The response is modeled by first projecting the stimulus onto each
% of the four filters, squaring their outputs, and taking the difference
% between the rightward and leftward detectors:
Eflts = [reshape(rightward_1,1,nX*nT)', reshape(rightward_2,1,nX*nT)'];
Sflts = [reshape(leftward_1,1,nX*nT)', reshape(leftward_2,1,nX*nT)'];
genpot = sum((allStim*Eflts)'.^2)' - sum((allStim*Sflts)'.^2)';

% The signal is now in units of firing rate.  We simulate spikes with an
% approximately Poisson process:
nspks = (rand(size(genpot,1),1) < genpot);
spikes = [zeros(1,nT-1), nspks(1:end-2)'];
sum(spikes)
stim = stim(1:length(spikes),:);

%% SPIKE-TRIGGERED COVARIANCE -------------------------------------------
%
% Above, we simulated a white noise experiment.  Now we will recover a
% model of our simulated neuron using spike-triggered covariance.

% We can't use the first few spikes because they don't have a complete
% stimulus history:
spikes(1:(nT-1)) = 0; 
nspikes = sum(spikes); 
spikeInd=(find(spikes>0.5));
wts = spikes(spikeInd);

% We begin by isolating the stimuli that produced spikes
spikeStim = zeros(length(spikeInd), nX*nT);
for n = 1:length(spikeInd)
  ind = spikeInd(n);
  spikeStim(n,:) = allStim(ind-nT+1,:);
end
 
% Now we compute the normalized (unit vector) STA and project it
% out of the spike-triggered stimuli.  Normally, one would compute covariance
% about a mean by subtracting the mean out (to recenter the distribution).  
% The filters recovered by STC will all be orthogonal;  projecting out the 
% STA ensures that the STA is orthogonal to these as well.  In this case,
% the STA will be flat (we are working with the phase-insensitive
% motion-energy model). But if we did have a residual STA, this would be
% important
figure;
sta = wts*spikeStim ./ nspikes;
staN = sta(:)/sqrt(sum(sta.^2)); % normalized STA
stim2 = spikeStim - spikeStim*staN*staN';  %project out the STA

% Now, we compute the covariance of stimuli preceding spikes:
covv = innerProd(stim2)./ (length(spikes)-1);
subplot(2,2,1);
% The diagonal of this matrix describes the variance, which swamps the
% covariance between dimensions. To visualize the covariance matrix, we
% will set the diagonal to zero so we can see the covariance structure:
imagesc(covv.*(1-eye(size(covv,1))));colormap(gray);
title('Covariance matrix')

% The dimensionality of our space (D) is defined as the number of bars (nX) 
% multiplied by the number of time bins (nT). The stimulus space can be 
% envisioned as a D dimensional space and each nX*nT stimulus block can be 
% envisioned as a vector in that space.  An "axis" in this space
% corresponds to a particular instantiation of the stimulus on once side of
% the origin, the inverse of that stimulus on the other side of the
% origin, and all points in between. The covariance matrix describes the 
% dependency of spiking on the covariation between every pair of dimensions 
% (e.g. how spiking depends on the covariation between bar 2 at time 2 and bar 
% 2 at time 3). A principal components analysis (PCA) applied to this matrix 
% returns D orthogonal vectors (the eigenvectors) and the variance of the 
% distribution spike-triggered stimulus distribution along each vectors, the 
% eigenvalues.   
[evect,eval] = sortedEig(covv);  eval = diag(eval);

% The stimuli that produced spikes form a cloud of points in this high
% dimensional space. PCA fits a D dimensional hyperellipse to this cloud.
% Given no correlation between the stimulus and the spikes and an infinite 
% amount of data, the vectors will be assigned randomly and the variance 
% along each axis will be the variance of the stimulus (PCA will fit a
% hypersphere to the spike-triggered stimuli).  Given no correlation
% between the stimulus and the spikes and a finite amount of data, some axes
% will have a variance slightly higher than the data and some lower, just by 
% chance.  We are interested in vectors along which the variance is significantly
% different than these chance correlations.  Normally, we would use bootstrap
% hypothesis testing to confirm which filters have a variance significantly
% different than the chance distribution.  But for our purposes here, a
% thumb-rule will suffice: when the filters are plotted in rank order,
% those that "pop-off" of the distribution are significant.  

% Plot the eigenvalues
subplot(2,2,2)
mxx = max(max(eval(1:end-1))).*1.3;
plot(eval(1:end-1),'.')
axis([0 length(evect) 0 mxx])
xlabel('rank number')
ylabel('eigenvalue (variance)')
title('PCA of the covariance matrix')

% Now, let's plot the filters.  First, plot the STA.  The STA is
% flat (approximately zero valued).
subplot(2,5,6)
mx = max([max(max(abs(staN))), max(max(abs(evect)))]);
imagesc(reshape(staN,nT,nX),[-mx mx]);
axis('off');title('STA'); colormap(gray);

% Plot the STC filters.  Remember, only those filters whose eigenvalues
% are significantly higher or lower than chance ("pop-off" the
% distribtuion) are considered significant.  If you look closely, you will
% see that there are two filters with high variance and two with low
% variance that meet this criteria:
filters = evect(:,[1 2 nT*nX-2 nT*nX-1]);
for n=1:size(filters,2)
    subplot(2,5,6+n)
    imagesc(reshape(filters(:,n),nT,nX),[-mx mx]);
    title(['STC Filter ',num2str(n)]);
    axis('off')
end;
colormap(gray)

% Compare these filters to the original rightward and leftward filters we
% started with.  Note that the filters may be inversed in sign relative to
% the original filters.  In our simulated neuron, the filter outputs were
% squared, thus it is irrelevant whether the filter or its inverse are returned by
% PCA.  Note also that PCA forces the filters to be orthogonal and consequently the
% filters need not be exactly those that we started with.  But the filters recovered 
% by STC span the same linear subspace as the filters for the simulated neuron.  
% In other words, the set of filters will produce a model that given the same 
% input (stimulus), will produce the same output.

% RECOVERING THE NONLINEARITY --------------------------------------------
%
% In our original model, the rightward and leftward detectors were each
% combined via a sum of squares and the difference between the the
% rightward and leftward signals taken to produce a firing rate.  Thus 
% far, we have recovered the set of linear filters.  Now we want to recover 
% the nonlinear function that describes the combination of filter outputs.  
% Given that we recovered 4 filters from STC, this nonlinear function is
% described by a joint 4-dimensional probability distribution.  Unfortunately, 
% we don't have enough data here to reconstruct this function completely. But
% we can look at slices through it.  Here, we reconstruct firing rate as a 
% function of the output of single filters by taken the ratio of the spiking 
% stimuli and all stimuli projected onto the filters.

figure
for n=1:4
    allProj = allStim*filters(:,n);
    spkProj = spikeStim*filters(:,n);

    mx = max(abs(allProj))+0.01;
    binedg = -mx:mx/19:mx;

    [rawhist] = histc(allProj, binedg);
    [spkhist] = histc(spkProj, binedg);

    NL(n,:) = spkhist(1:end-1)./rawhist(1:end-1);
    bins(n,:) = binedg(1:end-1)+diff(binedg(1:2))/2;
    subplot(4,1,n)
    plot(bins(n,:), NL(n,:));
    hold on;
    plot(bins(n,:), NL(n,:),'o');
    if (n<3)
    axis([-1 1 0 1])
    text(-0.2, 0.8,['STC filter ',num2str(n)]);
    else
    axis([-1 1 0 0.3])
    text(-0.2, 0.25,['STC filter ',num2str(n)]);
    end;
    ylabel('spikes/frame')
    
end;
    xlabel('filter output')
    
  % As you can see, firing rate increases as a function of the squared output of
  % STC filters 1 and 2 and decreases as a function of the squared output
  % of filters 3 and 4.  Thus, the first two filters have an excitatory
  % influence on the simulated neuron's response and the last two are
  % suppressive.  We now have a complete model of our simulated neuron that
  % describes how a stimulus is converted into a firing rate response.

% END