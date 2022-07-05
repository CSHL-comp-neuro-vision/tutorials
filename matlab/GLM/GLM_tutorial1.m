% GLM_tutorial1.m
%
% Author: Jonathan Pillow (June 18 2010)


addpath('code_GLM_v1_Feb2010\GLMcode');
addpath('code_GLM_v1_Feb2010\tools_misc');

addpath('code_GLM_v1_Feb2010\tools_mexcode');

%%
% This script provides a simple introduction to simulating and fitting a
% GLM point process model for neural spike trains.  For more complex
% examples, see the additional scripts in "testcripts/".
%
% If publishing scientific work making use of this code, please cite:
%
% Pillow JW, Shlens J, Paninski L, Sher A, Litke AM, Chichilnisky EJ,
% Simoncelli EP. (2008) Spatio-temporal correlations and visual signaling
% in a complete neuronal population. Nature 454: 995-999 
%
% The full GLM code for matlab can be downloaded from:
% http://pillowlab.cps.utexas.edu/code_GLM.html


% Basic summary of code in this script:
%   1. Set up model params
%      - 1D (temporal) stimulus filter
%      - post-spike filter with refractory and bursty effects
%      - exponential nonlinearity
%   2. Simulated model responses to white noise stimuli
%   3. Generate training data for fitting (Simulate response to a long stimulus)
%   4. Fit GLM params via maximum-likelihood (requires optimization toolbox)
%   5. Show comparison of true filter, STA and GLM maximum-likelihood estimates

%% ==== 1.Set up the GLM model (and other params for simulation) ====

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100; 

DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 20;    % Number of time bins in stimulus filter
ttk = [-nkt+1:0]';  % time relative to spike of stim filter taps
ggsim = makeSimStruct_GLM(nkt,DTsim); % Create GLM structure with default params


%% The above code creates a generic"GLM" structure with default parameters,
% but let's set some of them by hand to get a bit more intuition for what
% they are:

% 1a. Make stimulus filter
ggsim.k = normpdf(ttk, -3, 1.25) - .8*normpdf(ttk, -nkt/3, 4);  

clf; % Plot it
plot(ttk, ggsim.k);
title('stimulus filter');
ylabel('filter weight');
xlabel('time before spike (frames)');


%% 1b. Make post-spike filter
% Make basis for (h) 
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [0 3];  % Peak location for first and last vectors
ihbasprs.b = .5;  % How nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period (optional param)
[iht,ihbasOrthog,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);

% store these in the simulation structure
ggsim.ihbasprs = ihbasprs;  
ggsim.iht = iht;
%  (Note: ihbasOrthog is the orthogonal version of 'ihbasis', which makes
%  fitting slightly faster, but we'll worry about that later).

% Make post-spike filter h using this basis
ihweights = [-10 1 1 -1 -.5]';  % weights
ih = ihbasis*ihweights;  % post-spike filter
ggsim.ih = ih;

% ---- Plot them ----------
subplot(311);
plot(ggsim.iht, ihbasis);
title('basis for post-spike filter h');
subplot(312);
plot(ggsim.iht, ggsim.iht*0, 'k--', ggsim.iht, ggsim.ih);
title('post-spike filter h');
set(gca, 'ylim',[min(ggsim.ih)*1.1 max(ggsim.ih)*1.5]);
subplot(3,1,3); % --------
plot(ggsim.iht, exp(ggsim.iht*0), 'k--', ggsim.iht, exp(ggsim.ih));
title('exponentiated post-spike filter h');
xlabel('time after spike (frames)');
ylabel('gain');

% 1c. Finally, the additive "dc" input to the neuron
ggsim.dc = 3;  % (gives positive bias: baseline rate will be exp(ggsim.dc))

% 1d. Nonlinearity: keep to be exponential
ggsim.nlfun = @exp;


%% ==== 2. Make GWN stimulus & simulate the glm model response. =========

slen = 100; % Stimulus length (frames) 
swid = 1;  % Stimulus width  (pixels).  Must match # pixels in stim filter
Stim = randn(slen,swid);  % Gaussian white noise stimulus
[tsp, Itot,Ispk] = simGLM(ggsim, Stim);  % Simulate GLM response

% ==== Make Figure: repeat responses ========
tt = [DTsim:DTsim:slen]';
subplot(221); %------------------------
plot(1:slen, Stim, 'k', 'linewidth', 2); 
title('GWN stimulus');
xlabel('time');
axis tight;

subplot(222); %------------------------
plot(tt, Itot-Ispk, tt, Ispk, 'r');
title('stim (blue) and spike-hist (red) filter outputs'); axis tight;

subplot(223); %------------------------
if ~isempty(tsp)
    plot(tt, Itot, tsp, max(Itot)*ones(size(tsp)), 'r.');
else
    plot(tt,Itot);
end
title('net filter output (red) and spike times (dots)');
axis tight;
xlabel('time (frames)');


% -----Run a few repeat simulations ------
nrpts = 5;        % number of repeats to draw
subplot(224);
for j = 1:nrpts;
    [tsp1,Itot1] = simGLM(ggsim,Stim);
    plot(tt,Itot1, 'color', rand(3,1)); hold on;
end
axis tight; hold off;
title('repeated responses to same stimulus');
xlabel('time (frames)');

%% 2b.  The above is highly optimized simulation code: let's write our own
% simple version right here:

% Step 1: filter the stimulus with the stimulus filter:
Sfilt = sameconv(Stim,ggsim.k);  
subplot(211);
plot(Stim); title('raw stimulus');
subplot(212);
plot(Sfilt); title('filtered stimulus');

% Step 2: add dc term
Sfilt = Sfilt + ggsim.dc;

% Step 3: re-bin time into very small bins (so we can capture the fine
% precision of spike trains);
ttlow = (1:slen)';  % original time indices (in stimulus frames);
tthi = (.5+DTsim:DTsim:slen+.5)';  % fine timescale bins;

% Step 4: linearly interpolate the stimulus to these fine time bins
StimHi = interp1(ttlow, Sfilt, tthi, 'linear', 'extrap');  

% Step 5: finally, we need a loop to generate spikes and add in the
% post-spike filter contribution every time there's a spike
rlen = length(StimHi); % total length of high-frequency time bins
hlen = length(ggsim.ih);% length of post-spike filter
II = StimHi;  % This stores the total filter output (stim + spike filter)
sptimes = [];

% Loop: apply nonlinearity to net filter output and draw spikes
for j = 1:rlen
    % Draw a bernoulli random variable to determine whether to spike or not
    if rand < ggsim.nlfun(II(j))*DTsim/RefreshRate
        % Spike! 
        sptimes = [sptimes; DTsim*j];
        % add h to subsequent time bins
        II(j+1:min(rlen,j+hlen)) = II(j+1:min(rlen,j+hlen)) + ...
            ggsim.ih(1:length(j+1:min(rlen,j+hlen)));
    end
end

% That's all!  Plot it:
subplot(211);
plot(tthi, StimHi, tthi, II);
title('Stim filter + net filter output');
axis tight;
subplot(212);
plot(tthi, II-StimHi, sptimes, 1, 'r.');
title('h-filter output and spikes');
axis tight;


%% 3. Make some (longer) training data for fitting %======================
slen = 2500; % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,swid))*2-1;  %  Run model on long, binary stimulus
[tsp,Itot,ispk] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp);

% Compute STA and use as initial guess for k
sta0 = simpleSTC(Stim,tsp,nkt);
sta = reshape(sta0,nkt,[]);

% -----------
% Make param object with "true" params;
ggTrue = makeFittingStruct_GLM(ggsim.k,DTsim,ggsim);
ggTrue.tsp = tsp;
ggTrue.tspi = 1;  % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)

% Check that conditional intensity calculated under the simulation code
% (above) is the same as the conditional intensity computed under the
% log-likelihood code 

[neglogliTrue, rrT,tt] = neglogli_GLM(ggTrue,Stim);  % compute negative log-likelihood
nsamps = 10000;
clf;
plot(tt(1:nsamps), rrT(1:nsamps), tt(1:nsamps), exp(Itot(1:nsamps)), 'r--');
 
% Finally, show that we could compute the log-likelihood from our
% conditional intensity rrT and spike times very simply:
spkInds = round(tsp./DTsim);
logli = sum(log(rrT(spkInds))) - sum(rrT)*DTsim/RefreshRate; % Compute log-li
[logli -neglogliTrue]  % these should agree

% -----------

%% 4. Do ML fitting of params with simulated data %=====================

%  Initialize params for fitting --------------
gg0 = makeFittingStruct_GLM(sta,DTsim);  % projects sta into basis for fitting k
gg0.tsp = tsp;  % Insert spikes into fitting struct
gg0.tspi = 1;   % First spike to use (you can ask it to ignore the first "n" spikes)
gg0.iht = ggsim.iht; % set post-spike filter times and post-spike filter basis
gg0.ihbas = ihbasis;
[logli0,rr0,tt] = neglogli_GLM(gg0,Stim); % Compute logli of initial params (if desired)

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100};
[gg, negloglival] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% 5. Plot results ======================================================
ttk = -nkt+1:0;
subplot(221);  % True filter  % ---------------
plot(ttk, ggsim.k, 'k', ttk, sta, ttk, gg.k, 'r');
title('Stim filters estimates');
legend('true', 'STA', 'ML', 'location', 'northwest');

subplot(223);
flts = normalizecols([ggsim.k, sta, gg.k]);
plot(ttk, flts(:,1),'k', ttk,flts(:,2), ttk, flts(:,3), 'r');
title('Normalized Stim filters');
xlabel('time before spike (frames)');

subplot(222); % ----------------------------------
plot(ggsim.iht, ggsim.ih, 'k', gg.iht, gg.ihbas*gg.ih,'r');
title('post-spike kernel');
legend('true', 'ML', 'location', 'southeast');
axis tight;

subplot(224); % ----------------------------------
plot(ggsim.iht, exp(ggsim.ih), 'k', gg.iht, exp(gg.ihbas*gg.ih),'r');
title('exponentiated post-spike kernel');
xlabel('time since spike (frames)');
ylabel('gain');
axis tight;

% Errors in STA and ML estimate (compute the angle between the true and
% estimated filters)
Estim_Error = [subspace(flts(:,1),flts(:,2)), subspace(flts(:,1),flts(:,3))] % In radians
