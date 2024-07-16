% tutorial3_regularization_linGauss.m
%
% This is an interactive tutorial designed to walk through
% regularization for a linear-Gaussian GLM, which allows for closed-form
% MAP parameter estimates.  The next tutorial ('tutorial4') will cover the
% same methods for the Poisson GLM (which requires numerical optimization).
% 
% We'll consider two simple regularization methods:
%
% 1. Ridge regression - corresponds to maximum a posteriori (MAP)
%                       estimation under an iid Gaussian prior on the
%                       filter coefficients. 
%
% 2. L2 smoothing prior - using to an iid Gaussian prior on the
%                         pairwise-differences of the filter(s).
%
% Data: from Uzzell & Chichilnisky 2004; see README file for details. 
%
% Last updated: Mar 10, 2020 (JW Pillow)

% Tutorial instructions: Execute each section below separately using
% cmd-enter. For detailed suggestions on how to interact with this
% tutorial, see header material in tutorial1_PoissonGLM.m

%% ====  1. Load the raw data ============

% ------------------------------------------------------------------------
% Be sure to unzip the data file data_RGCs.zip
% (http://pillowlab.princeton.edu/data/data_RGCs.zip) and place it in 
% this directory before running the tutorial.  
% ------------------------------------------------------------------------
% (Data from Uzzell & Chichilnisky 2004):
datdir = 'data_RGCs/';  % directory where stimulus lives
load([datdir, 'Stim']);    % stimulus (temporal binary white noise)
load([datdir,'stimtimes']); % stim frame times in seconds (if desired)
load([datdir, 'SpTimes']); % load spike times (in units of stim frames)
ncells = length(SpTimes);  % number of neurons (4 for this dataset).
% Neurons #1-2 are OFF, #3-4 are ON.
% -------------------------------------------------------------------------

% Pick a cell to work with
cellnum = 3; % (1-2 are OFF cells; 3-4 are ON cells).
tsp = SpTimes{cellnum};

% Compute some basic statistics on the stimulus
dtStim = (stimtimes(2)-stimtimes(1)); % time bin size for stimulus (s)

% See tutorial 1 for some code to visualize the raw data!

%% == 2. Upsample to get finer timescale representation of stim and spikes === 

% The need to regularize GLM parameter estimates is acute when we don't
% have enough data relative to the number of parameters we're trying to
% estimate, or when using correlated (eg naturalistic) stimuli, since the
% stimuli don't have enough power at all frequencies to estimate all
% frequency components of the filter. 
%
% The RGC dataset we've looked at so far requires only a temporal filter
% (as opposed to spatio-temporal filter for full spatiotemporal movie
% stimuli), so it doesn't have that many parameters to esimate. It also has
% binary white noise stimuli, which have equal energy at all frequencies.
% Regularization thus isn't an especially big deal for this data (which was
% part of our reason for selecting it). However, we can make it look
% correlated by considering it on a finer timescale than the frame rate of
% the monitor.  (Indeed, this will make it look highly correlated).

% For speed of our code and to illustrate the advantages of regularization,
% let's use only a reduced (5-minute) portion of the dataset:
nT = 120*60*1;  % # of time bins for 1 minute of data 
Stim = Stim(1:nT); % pare down stimulus
tsp = tsp(tsp<nT*dtStim); % pare down spikes

% Now upsample to finer temporal grid
upsampfactor = 10; % divide each time bin by this factor
dtStimhi = dtStim/upsampfactor; % use bins 100 time bins finer
ttgridhi = (dtStimhi/2:dtStimhi:nT*dtStim)'; % fine time grid for upsampled stim
Stimhi = interp1((1:nT)*dtStim,Stim,ttgridhi,'nearest','extrap');
nThi = nT*upsampfactor;  % length of upsampled stimulus

% Visualize the new (upsampled) raw data:
subplot(211);
iiplot = 1:(60*upsampfactor); % bins of stimulus to plot
ttplot = iiplot*dtStimhi; % time bins of stimulus
plot(ttplot,Stimhi(iiplot), 'linewidth', 2);  axis tight;
title('raw stimulus (fine time bins)');
ylabel('stim intensity');
% Should notice stimulus now constant for many bins in a row

% Bin the spike train and replot binned counts
sps = hist(tsp,ttgridhi)';
subplot(212);
stem(ttplot,sps(iiplot));
title('binned spike counts');
ylabel('spike count'); xlabel('time (s)');
set(gca,'xlim', ttplot([1 end]));
% Now maximum 1 spike per bin!

%%  3. Divide data into "training" and "test" sets for cross-validation

trainfrac = .8;  % fraction of data to use for training
ntrain = ceil(nThi*trainfrac);  % number of training samples
ntest = nThi-ntrain; % number of test samples
iitest = 1:ntest; % time indices for test
iitrain = ntest+1:nThi;   % time indices for training
stimtrain = Stimhi(iitrain,:); % training stimulus
stimtest = Stimhi(iitest,:); % test stimulus
spstrain = sps(iitrain,:);
spstest =  sps(iitest,:);

fprintf('Dividing data into training and test sets:\n');
fprintf('Training: %d samples (%d spikes) \n', ntrain, sum(spstrain));
fprintf('    Test: %d samples (%d spikes)\n', ntest, sum(spstest));

%% === 4. Fit the linear-Gaussian model using ML ====================

% Set the number of time bins of stimulus to use for predicting spikes
ntfilt = 20*upsampfactor;  % Try varying this, to see how performance changes!

% build the design matrix, training data
Xtrain = [ones(ntrain,1),....
    hankel([zeros(ntfilt-1,1);stimtrain(1:end-ntfilt+1)],...
    stimtrain(end-ntfilt+1:end))];

% Compute maximum likelihood filter estimate ("whitened STA");
filtML = (Xtrain'*Xtrain)\(Xtrain'*spstrain);
ttk = (-ntfilt+1:0)*dtStimhi;
plot(ttk,filtML(2:end));
xlabel('time before spike');

% Looks like garbage! If you reduce 'upsampfactor' to 2, it still looks ok,
% but beyond this the stimulus lacks support at high frequencies and the
% covariance (Xdsgn'*Xdsgn) becomes badly conditioned.

%% === 5. Ridge regression: linear-Gaussian model ======================

% Now let's regularize by adding a penalty on the sum of squared filter
% coefficients w(i) of the form:   
%       penalty(lambda) = lambda*(sum_i w(i).^2),
% where lambda is known as the "ridge" parameter.  This is also known as an
% "L2 penalty".  Minimizing error plus this penalty ("penalized least
% squares") is equivalent to computing the MAP estimate under an iid
% Gaussian prior on the filter coefficients.  
%
% The MAP estimate for the LG model parameters has a closed form, making it
% simple and fast to compute: w_hat = (X'*X + lambda*I)^{-1} * X^T*Y 
%
% The only remaining question is: how to set lambda?
% We'll show here how to do it for a grid of lambda values and use
% cross-validation (test error) to select which is best.  


% Set up grid of lambda values (ridge parameters)
lamvals = 2.^(0:15); % it's common to use a log-spaced set of values
nlam = length(lamvals);

% Build design matrix for test data
Xtest = [ones(ntest,1), hankel([zeros(ntfilt-1,1); stimtest(1:end-ntfilt+1)], ...
    stimtest(end-ntfilt+1:end))];

% Precompute some quantities (X'X and X'*y) for training and test data
XXtr = Xtrain'*Xtrain;
XYtr = Xtrain'*spstrain;  % spike-triggered average, training data
Imat = eye(ntfilt+1); % identity matrix of size of filter + const
Imat(1,1) = 0; % don't apply penalty to constant coeff

% Allocate space for train and test errors
msetrain = zeros(nlam,1);  % training error
msetest = zeros(nlam,1);   % test error
w_ridge = zeros(ntfilt+1,nlam); % filters for each lambda

% Now compute MAP estimate for each ridge parameter
clf; plot(ttk,ttk*0,'k-'); hold on;
for jj = 1:nlam
    
    % Compute ridge regression estimate
    w = (XXtr+lamvals(jj)*Imat) \ XYtr; 
    
    % Compute MSE
    msetrain(jj) = (mean((spstrain-Xtrain*w).^2)); % training error
    msetest(jj) = (mean((spstest-Xtest*w).^2)); % test error
    
    % store the filter
    w_ridge(:,jj) = w;
    
    % plot it
    plot(ttk,w(2:end));
    title(['ridge estimate: lambda = ', num2str(lamvals(jj))]);
    xlabel('time before spike (s)'); drawnow; pause(.2);
 
end
hold off;
% note that the esimate "shrinks" down as we increase lambda

%% Plot filter estimates and errors for ridge estimates

subplot(222);
plot(ttk,w_ridge(2:end,:)); axis tight;  
title('all ridge estimates');
subplot(221);
semilogx(lamvals,msetrain,'o-', 'linewidth', 2);
title('training error');
subplot(223); 
semilogx(lamvals,msetest,'-o', 'linewidth', 2);
xlabel('lambda');
title('test error');

% Notice that training error gets monotonically worse as we increase lambda
% However, test error has an dip at some optimal, intermediate value.

% Determine which lambda is best by selecting one with lowest test error 
[~,imin] = min(msetest);
filt_ridge= w_ridge(2:end,imin);
subplot(224);
plot(ttk,ttk*0, 'k--', ttk,filt_ridge,'linewidth', 2);
xlabel('time before spike (s)'); axis tight;
title('best ridge estimate');


%% === 6. L2 smoothing: linear-Gaussian model ===========================

% Now let's instead instead try putting a penalty (or prior) on the squared
% differences between filter coefficients. This penalize large jumps in the
% filter, encouraging smoothness.

% First, we need to make a matrix D such that w'*D*w computes the squared
% differences. We can do this fairly easy as follows:

% This matrix computes differences between adjacent coeffs
Dx1 = spdiags(ones(ntfilt,1)*[-1 1],0:1,ntfilt-1,ntfilt); 
Dx = Dx1'*Dx1; % computes squared diffs

subplot(121); 
imagesc(Dx1(1:20,1:20)); axis image; title('computes diffs');
subplot(122);
imagesc(Dx(1:20,1:20)); axis image; title('computes squared diffs');

% Let's just check to be sure
x = randn(ntfilt,1); % make random filter
fprintf('sum of squared diffs, direct way: %.2f\n', sum(diff(x).^2));
fprintf('sum of squared diffs, matrix way: %.2f\n', x'*Dx*x);

%% Select smoothing penalty by cross-validation 

% Now let's do exactly what we did with a grid of lambda values, only
% instead of a scaled identity matrix we'll be using a scaled version of
% the 'Dx' matrix.

% Set up grid of lambda values (ridge parameters)
lamvals = 2.^(5:22); % it's common to use a log-spaced set of values
nlam = length(lamvals);

% Embed Dx matrix in matrix with one extra row/column for constant coeff
D = blkdiag(0,Dx); 

% Allocate space for train and test errors
msetrain_sm = zeros(nlam,1);  % training error
msetest_sm = zeros(nlam,1);   % test error
w_smooth = zeros(ntfilt+1,nlam); % filters for each lambda

% Now compute MAP estimate for each ridge parameter
clf; plot(ttk,ttk*0,'k-'); hold on;
for jj = 1:nlam
    
    % Compute ridge regression estimate
    w = (XXtr+lamvals(jj)*D) \ XYtr; 
    
    % Compute MSE
    msetrain_sm(jj) = mean((spstrain-Xtrain*w).^2); % training error
    msetest_sm(jj) =  mean((spstest-Xtest*w).^2); % test error
    
    % store the filter
    w_smooth(:,jj) = w;
    
    % plot it
    plot(ttk,w(2:end)); 
    title(['smoothing estimate: lambda = ', num2str(lamvals(jj))]);
    xlabel('time before spike (s)'); drawnow; pause(.2);
 
end
hold off;
% Note that the esimate becomes smoother as lambda increases.

%% Plot filter estimates and errors for smoothing estimates

subplot(222);
plot(ttk,w_smooth(2:end,:)); axis tight;  
title('all smoothing estimates');
subplot(221);
semilogx(lamvals,msetrain_sm,'o-', 'linewidth', 2);
title('training error');
subplot(223); 
semilogx(lamvals,msetest_sm,'-o', 'linewidth', 2);
xlabel('lambda');
title('test error');

% Notice that training error gets monotonically worse as we increase lambda
% However, test error has an dip at some optimal, intermediate value.

% Determine which lambda is best by selecting one with lowest test error 
[~,imin] = min(msetest_sm);
filt_smooth= w_smooth(2:end,imin);
subplot(224);
h = plot(ttk,ttk*0, 'k--', ttk,filt_ridge,...
    ttk,filt_smooth,'linewidth', 2);
xlabel('time before spike (s)'); axis tight;
title('best smoothing estimate');
legend(h(2:3), 'ridge', 'L2 smoothing', 'location', 'northwest');
% clearly the "L2 smoothing" filter looks better by eye!

% Last, lets see which one actually achieved lower test error
fprintf('\nBest ridge test error:      %.5f\n', min(msetest));
fprintf('Best smoothing test error:  %.5f\n', min(msetest_sm));

% These differences are pretty puny for this dataset, but they are
% greater if we examine Poisson-GLM log-likelihood instead of the
% linear-Gaussian model considered here. (See tutorial4).

