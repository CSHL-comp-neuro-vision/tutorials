% tutorial1b_PoissonGLM.m
%
% This is a tutorial illustrating the fitting of a linear-Gaussian GLM
% (also known as linear least-squares regression model) and a Poisson GLM
% (aka  "linear-nonlinear-Poisson" model) to retinal ganglion cell spike
% trains stimulated with binary temporal white noise. 
%
% NOTE: this is a variant of `tutorial1_PoissonGLM` which was altered just
% to remove stimulus filters.  So the model here uses spike-history only.
%
%
%
% DATASET: this tutorial is designed to run with retinal ganglion cell data
% from Uzzell & Chichilnisky 2004. Unfortunately, data is not provided with
% the public version of this repository because the dataset is not yet open
% access. If you would like to access the dataset, please write to me
% directly (pillow at princeton.edu).
%
% Last updated: July 18, 2022 (JW Pillow)


%% ====  1. Load the raw data ============

% ------------------------------------------------------------------------
% See README.md for info about gaining access to this dataset, or
% substitute your own data here. 
% ------------------------------------------------------------------------

% (Data from Uzzell & Chichilnisky 2004):
datdir = 'data_RGCs/';  % directory where stimulus lives
load([datdir, 'Stim']);    % stimulus (temporal binary white noise)
load([datdir,'stimtimes']); % stim frame times in seconds (if desired)
load([datdir, 'SpTimes']); % load spike times (in units of stim frames)

% Pick a cell to work with
cellnum = 1; % (1-2 are OFF cells; 3-4 are ON cells).
tsp = SpTimes{cellnum};
% -------------------------------------------------------------------------

% Compute some basic statistics on the data
dtStim = (stimtimes(2)-stimtimes(1)); % time bin size for stimulus (s)
RefreshRate = 1/dtStim; % Refresh rate of the monitor
nT = size(Stim,1); % number of time bins in stimulus
nsp = length(tsp); % number of spikes

% Print out some basic info
fprintf('--------------------------\n');
fprintf('Loaded RGC data: cell %d\n', cellnum);
fprintf('Number of stim frames: %d  (%.1f minutes)\n', nT, nT*dtStim/60);
fprintf('Time bin size: %.1f ms\n', dtStim*1000);
fprintf('Number of spikes: %d (mean rate=%.1f Hz)\n\n', nsp, nsp/nT*RefreshRate);

% Let's visualize some of the raw data
subplot(211);
iiplot = 1:120; % bins of stimulus to plot
ttplot = iiplot*dtStim; % time bins of stimulus
plot(ttplot,Stim(iiplot), 'linewidth', 2);  axis tight;
title('raw stimulus (full field flicker)');
ylabel('stim intensity');
subplot(212);
tspplot = tsp((tsp>=ttplot(1))&(tsp<ttplot(end)));
plot(tspplot, 1, 'ko', 'markerfacecolor', 'k');
set(gca,'xlim', ttplot([1 end]));
title('spike times'); xlabel('time (s)');

%% ==== 2. Bin the spike train ===== 
%
% For now we will assume we want to use the same time bin size as the time
% bins used for the stimulus. Later, though, we'll wish to vary this.

tbins = (.5:nT)*dtStim; % time bin centers for spike train binnning
ncells = length(SpTimes);  % number of neurons (4 for this dataset).

%sps = hist(tsp,tbins)';  % binned spike train

% For now we will assume we want to use the same time bin size as the time
% bins used for the stimulus. Later, though, we'll wish to vary this.
tbins = (.5:nT)*dtStim; % time bin centers for spike train binnning
sps = zeros(nT,ncells);
for jj = 1:ncells
    sps(:,jj) = hist(SpTimes{jj},tbins)';  % binned spike train
end
sps = sps(:,cellnum);

% Replot the responses we'll putting into our regression as counts
subplot(212);
stem(ttplot,sps(iiplot), 'k', 'linewidth', 2);
title('binned spike counts');
ylabel('spike count'); xlabel('time (s)');
set(gca,'xlim', ttplot([1 end]), 'ylim', [0 3.5]);


%% ==== 3. Build the design matrix: slow version ======
% This is a necessary step before we can fit the model: assemble a matrix
% that contains the relevant regressors for each time bin of the response,
% known as a design matrix.  Each row of this matrix contains the relevant
% stimulus chunk for predicting the spike count at a given time bin

% Set the number of time bins of stimulus to use for predicting spikes
ntfilt = 25;  % Try varying this, to see how performance changes!
nthist = ntfilt;


% Build spike-history design matrix

paddedSps = [zeros(nthist,1); sps(1:end-1,1)];
% SUPER important: note that this doesn't include the spike count for the
% bin we're predicting? The spike train is shifted by one bin (back in
% time) relative to the stimulus design matrix
Xdsgn = hankel(paddedSps(1:end-nthist+1), paddedSps(end-nthist+1:end));

% % Build design matrix without using a for loop
% paddedStim = [zeros(ntfilt-1,1); Stim]; % pad early bins of stimulus with zero
% Xdsgn = hankel(paddedStim(1:end-ntfilt+1), Stim(end-ntfilt+1:end));

% (You can check for you like that this gives the same matrix as the one
% created above!)

imagesc(-ntfilt+1:0, 1:50, Xdsgn(1:50,:));
xlabel('lags before spike time bin');
ylabel('time bin of response');
title('Design matrix');


%% === 4. Compute and visualize the spike-triggered average (STA) ====

% When the stimulus is Gaussian white noise, the STA provides an unbiased
% estimator for the filter in a GLM / LNP model (as long as the nonlinearity
% results in an STA whose expectation is not zero; feel free 
% to ignore this parenthetical remark if you're not interested in technical
% details. It just means that if the nonlinearity is symmetric, 
% eg. x^2, then this condition won't hold, and the STA won't be useful).
%
% In many cases it's useful to visualize the STA (even if your stimuli are
% not white noise), just because if we don't see any kind of structure then
% this may indicate that we have a problem (e.g., a mismatch between the
% design matrix and binned spike counts.

% It's extremely easy to compute the STA now that we have the design matrix
sta = (Xdsgn'*sps)/nsp;

% Plot it
ttk = (-ntfilt+1:0)*dtStim; % time bins for STA (in seconds)
plot(ttk,ttk*0, 'k--', ttk, sta, 'o-', 'linewidth', 2); axis tight;
title('STA'); xlabel('time before spike (s)');

% If you're still using cell #1, this should look like a biphasic filter
% with a negative lobe just prior to the spike time.

% (By contrast, if this looks like garbage then it's a good chance we did
% something wrong!)

%% 4b. whitened STA (ML fit to filter for linear-Gaussian GLM)

% If the stimuli are non-white, then the STA is generally a biased
% estimator for the linear filter. In this case we may wish to compute the
% "whitened" STA, which is also the maximum-likelihood estimator for the filter of a 
% GLM with "identity" nonlinearity and Gaussian noise (also known as
% least-squares regression). 
%
% If the stimuli have correlations this ML estimate may look like garbage
% (more on this later when we come to "regularization").  But for this
% dataset the stimuli are white, so we don't (in general) expect a big
% difference from the STA.  (This is because the whitening matrix
% (Xdsng'*Xdsgn)^{-1} is close to a scaled version of the identity.

% whitened STA
wsta = (Xdsgn'*Xdsgn)\sta*nsp;
% or equivalently inv(Xdsgn'*Xdsgn)*(Xdsgn'*sps)
% this is just the least-squares regression formula!

% Let's plot them both (rescaled as unit vectors so we can see differences
% in their shape).
h = plot(ttk,ttk*0, 'k--', ttk, sta./norm(sta), 'o-',...
    ttk, wsta./norm(wsta), 'o-', 'linewidth', 2); axis tight;
legend(h(2:3), 'STA', 'wSTA', 'location', 'northwest');
title('STA and whitened STA'); xlabel('time before spike (s)');

%% 4c. Predicting spikes with a linear-Gaussian GLM

% The whitened STA can actually be used to predict spikes because it
% corresponds to a proper estimate of the model parameters (i.e., for a
% Gaussian GLM). Let's inspect this prediction

sppred_lgGLM = Xdsgn*wsta;  % predicted spikes from linear-Gaussian GLM

% Let's see how good this "prediction" is
% (Prediction in quotes because we are (for now) looking at the performance
% on training data, not test data... so it isn't really a prediction!)

% Plot real spike train and prediction
stem(ttplot,sps(iiplot)); hold on;
plot(ttplot,sppred_lgGLM(iiplot),'linewidth',2); hold off;
title('linear-Gaussian GLM: spike count prediction');
ylabel('spike count'); xlabel('time (s)');
set(gca,'xlim', ttplot([1 end]));
legend('spike count', 'lgGLM');

%% 4d. Fitting and predicting with a linear-Gaussian-GLM with offset

% Oops, one thing we forgot above was to include an offset or "constant"
% term in the design matrix, which will allow our prediction to have
% non-zero mean (since the stimulus here was normalized to have zero mean).

% Updated design matrix
Xdsgn2 = [ones(nT,1), Xdsgn]; % just add a column of ones

% Compute whitened STA
MLwts = (Xdsgn2'*Xdsgn2)\(Xdsgn2'*sps); % this is just the LS regression formula
const = MLwts(1); % the additive constant
wsta2 = MLwts(2:end); % the linear filter part

% Now redo prediction (with offset)
sppred_lgGLM2 = const + Xdsgn*wsta2;

% Plot this stuff
stem(ttplot,sps(iiplot)); hold on;
h=plot(ttplot,sppred_lgGLM(iiplot),ttplot,sppred_lgGLM2(iiplot)); 
set(h, 'linewidth', 2);  hold off;
title('linear-Gaussian GLM: spike count prediction');
ylabel('spike count'); xlabel('time (s)');
set(gca,'xlim', ttplot([1 end]));
legend('spike count', 'lgGLM', 'lgGLM w/ offset');

% Let's report the relevant training error (squared prediction error on
% training data) so far just to see how we're doing:
mse1 = mean((sps-sppred_lgGLM).^2); % mean squared error, GLM no offset
mse2 = mean((sps-sppred_lgGLM2).^2);% mean squared error, with offset
rss = mean((sps-mean(sps)).^2); % squared error of spike train
fprintf('Training perf (R^2): lin-gauss GLM, no offset: %.2f\n',1-mse1/rss);
fprintf('Training perf (R^2): lin-gauss GLM, w/ offset: %.2f\n',1-mse2/rss);

%% ======  5. Poisson GLM ====================

% Let's finally move on to the LNP / Poisson GLM!

% This is super-easy if we rely on built-in GLM fitting code
pGLMwts = glmfit(Xdsgn,sps,'poisson', 'constant', 'on');
pGLMconst = pGLMwts(1); % constant ("dc term"); 
pGLMfilt = pGLMwts(2:end); % stimulus filter

% The 'glmfit' function will fit a GLM for us. Here we have specified that
% we want the noise model to be Poisson. The default setting for the link
% function (the inverse of the nonlinearity) is 'log', so default
% nonlinearity is 'exp'). The last argument tells glmfit to include a
% constant offset parameter in the model, so we can use the 'Xdsgn' matrix
% instead of the 'Xdsgn2' matrix that contains 1's in the first row.
% 
% The above call is thus equivalent to having done:
% > pGLMwts = glmfit(Xdsgn2,sps,'poisson', 'link', 'log','constant','off');

% Compute predicted spike rate on training data
ratepred_pGLM = exp(pGLMconst + Xdsgn*pGLMfilt);
%  (equivalent to if we had just written exp(Xdsgn2*pGLMwts))/dtStim;

%%  5b. Make plots showing and spike rate predictions

subplot(211);
h = plot(ttk,ttk*0, 'k--', ttk, wsta2./norm(wsta2), 'o-',...
    ttk, pGLMfilt./norm(pGLMfilt), 'o-', 'linewidth', 2); axis tight;
legend(h(2:3), 'lin-gauss GLM', 'poisson GLM', 'location', 'northwest');
title('(normalized) linear-Gaussian and Poisson GLM filter estimates'); 
xlabel('time before spike (s)');

subplot(212);
stem(ttplot,sps(iiplot)); hold on;
h = plot(ttplot,sppred_lgGLM2(iiplot),ttplot,ratepred_pGLM(iiplot)); 
set(h, 'linewidth', 2);  hold off;
title('spike rate predictions');
ylabel('spikes / bin'); xlabel('time (s)');
set(gca,'xlim', ttplot([1 end]));
legend('spike count', 'lin-gauss GLM', 'exp-poisson GLM');

% Note the rate prediction here is in units of spikes/bin. If we wanted
% spikes/sec, we could divide it by bin size dtStim.


%% 6. Non-parametric estimate of the nonlinearity

% The above fitting code assumes a GLM with an exponential nonlinearity
% (i.e., governing the mapping from filter output to instantaneous spike
% rate). We might wish to examine the adequacy of that assumption and make
% a "nonparametric" estimate of the nonlinearity using a more flexible
% class of functions.

% Let's use the family of piece-wise constant functions, which results in a
% very simple estimation procedure:
% 1. Bin the filter outputs
% 2. In each bin, compute the fraction of stimuli elicted spikes

% number of bins for parametrizing the nonlinearity f. (Try varying this!) 
nfbins = 30; 

% compute filtered stimulus
rawfilteroutput = pGLMconst + Xdsgn*pGLMfilt;

% bin filter output and get bin index for each filtered stimulus
[cts,binedges,binID] = histcounts(rawfilteroutput,nfbins); 
fx = binedges(1:end-1)+diff(binedges(1:2))/2; % use bin centers for x positions

% now compute mean spike count in each bin
fy = zeros(nfbins,1); % y values for nonlinearity
for jj = 1:nfbins
    fy(jj) = mean(sps(binID==jj));
end
fy = fy/dtStim; % divide by bin size to get units of sp/s;

% Now let's embed this in a function we can evaluate at any point
fnlin = @(x)(interp1(fx,fy,x,'nearest','extrap'));

% Make plots
subplot(211); % Plot exponential and nonparametric nonlinearity estimate
bar(fx,cts, 'hist');
ylabel('count'); title('histogram of filter outputs');

subplot(212);
xx = binedges(1):.01:binedges(end);
plot(xx,exp(xx)/dtStim,xx,fnlin(xx),'linewidth', 2);
xlabel('filter output');
ylabel('rate (sp/s)');
legend('exponential f', 'nonparametric f', 'location', 'northwest');
title('nonlinearity');

% What do you think of the exponential fit? Does this look like a good
% approximation to the nonparametric estimate of the nonlinearity?  Can you
% propose a better parametric nonlinearity to use instead?  

%
% Advanced exercise: write your own log-likelihood function that allows you
% to jointly optimize log-likelihood for the filter parameters and
% nonlinearity.  (By contrast, here we have optimized filter params under
% exponential nonlinearity and THEN fit the nonlinearity using these fixed
% filter parameters).  We could, for example, iteratively climb the
% log-likelihood as a function of filter params and nonlinearity params;
% this is a method known as "coordinate ascent").

%% 7. Quantifying performance: log-likelihood

% Lastly, compute log-likelihood for the Poisson GLMs we've used so far and
% compare performance.

% LOG-LIKELIHOOD (this is what glmfit maximizes when fitting the GLM):
% --------------
% Let s be the spike count in a bin and r is the predicted spike rate
% (known as "conditional intensity") in units of spikes/bin, then we have:   
%
%        Poisson likelihood:      P(s|r) = r^s/s! exp(-r)  
%     giving log-likelihood:  log P(s|r) =  s log r - r   
%
% (where we have ignored the -log s! term because it is independent of the
% parameters). The total log-likelihood is the summed log-likelihood over
% time bins in the experiment.

% 1. for GLM with exponential nonlinearity
ratepred_pGLM = exp(pGLMconst + Xdsgn*pGLMfilt); % rate under exp nonlinearity
LL_expGLM = sps'*log(ratepred_pGLM) - sum(ratepred_pGLM);

% 2. for GLM with non-parametric nonlinearity
ratepred_pGLMnp = dtStim*fnlin(pGLMconst + Xdsgn*pGLMfilt); % rate under nonpar nonlinearity
LL_npGLM = sps(sps>0)'*log(ratepred_pGLMnp(sps>0)) - sum(ratepred_pGLMnp);

% Now compute the rate under "homogeneous" Poisson model that assumes a
% constant firing rate with the correct mean spike count.
ratepred_const = nsp/nT;  % mean number of spikes / bin
LL0 = nsp*log(ratepred_const) - nT*ratepred_const;

% Single-spike information:
% ------------------------
% The difference of the loglikelihood and homogeneous-Poisson
% loglikelihood, normalized by the number of spikes, gives us an intuitive
% way to compare log-likelihoods in units of bits / spike.  This is a
% quantity known as the (empirical) single-spike information.
% [See Brenner et al, "Synergy in a Neural Code", Neural Comp 2000].
% You can think of this as the number of bits (number of yes/no questions
% that we can answer) about the times of spikes when we know the spike rate
% output by the model, compared to when we only know the (constant) mean
% spike rate. 

SSinfo_expGLM = (LL_expGLM - LL0)/nsp/log(2);
SSinfo_npGLM = (LL_npGLM - LL0)/nsp/log(2);
% (if we don't divide by log 2 we get it in nats)

fprintf('\n empirical single-spike information:\n ---------------------- \n');
fprintf('exp-GLM: %.2f bits/sp\n',SSinfo_expGLM);
fprintf(' np-GLM: %.2f bits/sp\n',SSinfo_npGLM);

% Let's plot the rate predictions for the two models 
% --------------------------------------------------
subplot(111);
stem(ttplot,sps(iiplot)); hold on;
plot(ttplot,ratepred_pGLM(iiplot),ttplot,ratepred_pGLMnp(iiplot),'linewidth',2); 
hold off; title('rate predictions');
ylabel('spikes / bin'); xlabel('time (s)');
set(gca,'xlim', ttplot([1 end]));
legend('spike count', 'exp-GLM', 'np-GLM');


%% 8. Quantifying performance: AIC

% Akaike information criterion (AIC) is a method for model comparison that
% uses the maximum likelihood, penalized by the number of parameters.
% (This allows us to compensate for the fact that models with more
% parameters can in general achieve higher log-likelihood. AIC determines
% how big this tradeoff should be in terms of the quantity:
%        AIC = - 2*log-likelihood + 2 * number-of-parameters
% The model with lower AIC is 
% their likelihood (at the ML estimate), penalized by the number of parameters  

AIC_expGLM = -2*LL_expGLM + 2*(1+ntfilt); 
AIC_npGLM = -2*LL_npGLM + 2*(1+ntfilt+nfbins);

fprintf('\n AIC comparison:\n ---------------------- \n');
fprintf('exp-GLM: %.1f\n',AIC_expGLM);
fprintf(' np-GLM: %.1f\n',AIC_npGLM);
fprintf('\nAIC diff (exp-np)= %.2f\n',AIC_expGLM-AIC_npGLM);
if AIC_expGLM < AIC_npGLM
    fprintf('AIC supports exponential-nonlinearity!\n');
else
    fprintf('AIC supports nonparametric nonlinearity!\n');
    % (despite its greater number of parameters)
end

% Caveat: technically the AIC should be using the maximum of the likelihood
% for a given model.  Here we actually have an underestimate of the
% log-likelihood for the non-parameteric nonlinearity GLM because
% because we left the filter parameters unchanged from the exponential-GLM.
% So a proper AIC comparison (i.e., if we'd achieved a true ML fit) would
% favor the non-parametric nonlinearity GLM even more!

% Exercise: go back and increase 'nfbins', the number of parameters (bins)
% governing the nonparametric nonlinearity. If you increase it enough, you
% should be able to flip the outcome so exponential nonlinearity wins.

% (Note: in the third tutorial we'll use cross-validation to properly
% evaluate the goodness of the fit of the models, e.g., allowing us to
% decide how many bins of stimulus history or how many bins to use for the
% non-parametric nonlinearity, or how to set regularization
% hyperparameters. The basic idea is to split data into training and test
% sets.  Fit the parameters on the training set, and compare models by
% evaluating log-likelihood on test set.)

%% 9. Simulating the GLM / making a raster plot

% Lastly, let's simulate the response of the GLM to a repeated stimulus and
% make raster plots 

iiplot = 1:60; % time bins of stimulus to use
ttplot = iiplot*dtStim; % time indices for these stimuli
StimRpt = Stim(iiplot); % repeat stimulus 
nrpts = 50;  % number of repeats
frate = exp(pGLMconst+Xdsgn(iiplot,:)*pGLMfilt);% firing rate in each bin

% Or uncomment this line to use the non-parametric nonlinearity instead:
%frate = dtStim*fnlin(pGLMconst+Xdsgn(iiplot,:)*pGLMfilt);% firing rate in each bin

% First, plot stimulus and true spikes
subplot(611);
plot(ttplot,Stim(iiplot), 'linewidth', 2);  axis tight;
title('raw stimulus (full field flicker)');
ylabel('stim intensity'); set(gca,'xticklabel', {});
subplot(612);
tspplot = tsp((tsp>=ttplot(1))&(tsp<ttplot(end)));
plot(tspplot, 1, 'ko', 'markerfacecolor', 'k');
set(gca,'xlim', ttplot([1 end]));
title('true spike times');
set(gca,'xticklabel', {});

% Simulate spikes using draws from a Bernoulli (coin flipping) process
spcounts = poissrnd(repmat(frate',nrpts,1)); % sample spike counts for each time bin
subplot(6,1,3:6);
imagesc(ttplot,1:nrpts, spcounts);
ylabel('repeat #');
xlabel('time (s)');
title('GLM spike trains');

%% Optional: redo using finer time bins, so we get maximum 1 spike per bin

upsampfactor = 100; % divide each time bin by this factor
dt_fine = dtStim/upsampfactor; % use bins 100 time bins finer
tt_fine = dt_fine/2:dt_fine:ttplot(end);

% Compute the fine-time-bin firing rate (which must be scaled down by bin width)
frate_fine = interp1(ttplot,frate,tt_fine,'nearest','extrap')'/upsampfactor;

% now draw fine-timescale spike train
spcounts_fine = poissrnd(repmat(frate_fine',nrpts,1)); % sample spike counts for each time bin

% Make plot
subplot(6,1,3:6);
imagesc(ttplot,1:nrpts, spcounts_fine);
ylabel('repeat #');
xlabel('time (s)');
title('GLM spike trains');


%% Suggested Exercises (advanced)
% -------------------------------
%
% 1) Go back and try it out for the other three neurons!
% (Go to block 1 and change the variable 'cellnum' to 1, 2, or 3.)
% 
% 2) Write your own code to do maximum likelihood estimation of the filter
% and the nonlinearity.  Your function should take in the parameters for
% the filter and the nonlinearity, and compute the Poisson log-likelihood
% function, the log Probability of the spike responses given the stimuli
% and the parameters.  A nice way to parametrize the nonlinearity is with a
% linear combination of basis functions, e.g.
%      f(x) = sum_i  w_i * f_i(x)
% where f_i(x) is the i'th basis function and w_i is the weight on that
% basis function.  You can choose the f_i to be Gaussian bumps or sigmoids,
% i.e. f_i(x) = 1./(1+exp(-x - c_i)) where c_i is the shift for the i'th
% basis function.
%
% Another alternative (that will prevent negative firing rates) is to
% parameterize the log-firing rate with a linear combination of basis
% functions, e.g.
%  log(f(x)) = sum_i  w_i * f_i(x)
%        meaning that
%  f(x) = exp(sum_i  w_i * f_i(x))
%  Now your weights can be negative or positive without fear of generating
%  negative values (which will cause your negative log-likelihood function
%  to give nans or -infs.  
%  
%  Write a function that takes in k and weight vector w and computes the
%  Poisson log-likelihood.  Hand that function off to fminunc and compare
%  the accuracy of the fits you get to the model with fixed exponential
%  nonlinearity.

