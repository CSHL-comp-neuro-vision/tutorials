%%DiffusionProcessTutorial.m
%
% This tutorial simulates the basic diffusion (random walk) model that
% predicts accuracy and RTs for a 2AFC decision.
%
% It requires adding the folder 'DiffusionProcess' to the path
%
% Written by G.M. Boynton, Summer 2008 For a great reference, see: Palmer,
% Huk & Shadlen (2008) Journal of Vision
%
% http://www.journalofvision.org/5/5/1/

addpath('dependencies_DiffusionProcess');

%% What is a 'diffusion process'?

% The basic idea is that decisions are made by accumulating information
% over time.  Specifically, evidence is represented as a variable that
% increments or decrements over time based on incoming information. I can't
% summarize it better than Palmer, Huk and Shadlen:

%  "The internal representation of the relevant stimulus is assumed to be
%  noisy and to vary over time. Each decision is based on repeated sampling
%  of this representation and comparing some function of these samples to a
%  criterion. For example, suppose samples of the noisy signal are taken at
%  discrete times and are added together to represent the evidence
%  accumulated over time. This accumulated evidence is compared to an upper
%  and lower bound. Upon reaching one of these bounds, the appropriate
%  response is initiated. If such a random walk model is modified by
%  reducing the time steps and evidence increments to infinitesimals, then
%  the model in continuous time is called a diffusion model (Ratcliff,
%  1978; Smith, 1990). For this model, the accumulated evidence has a
%  Gaussian distribution, which makes it a natural generalization of the
%  Gaussian version of signal detection theory (Ratcliff, 1980)."

%  (Diffusion processes are also known as the Wiener process or Brownian
%  motion)

% Intuitively, you can imagine how this simple model can predict error
% rates, and RTs for both correct and incorrect responses.  In the simplest
% case, analytical solutions exist for accuracy and the distribution of
% RTs.  

% This tutorial simulates a simple random walk model, calculating accuracy
% and histograms of RTs.  The results of the simulation is then compared to
% curves generated from the analytic solution.

%% Define variables for the random walk:

p.a = .1;  %upper bound (correct answer)
p.b = .1; %lower bound (incorrect answer)
p.u = .1; %drift rate (units/sec)
p.s = .1;  %standard deviation of drift (units/sec)


p.dt = .1; %step size for simulations (seconds)

%% Quick simulation

% First, let's do a quick simulation of the random walk model without
% worrying about keeping track of RTs.

nReps = 5; %number of staircases per simulation
nSteps = 30; %number of time steps

% Generate the entire matrix of step sizes in a single line. In this case,
% step sizes are pulled from a normal distribution with mean u*dt and
% standard deviation s*sqrt(dt).  Why sqrt(dt)?  Because variance adds
% linearly over time, so standard deviation adds by sqrt(dt).

dy = p.u*p.dt + p.s*sqrt(p.dt)*randn(nSteps,nReps);

% The random walk is a cumuluative sum of the steps (dy)
y = cumsum(dy);

% We need a time vector for plotting:
t = p.dt:p.dt:p.dt*nSteps;

% plot it.
figure(1)
clf
subplot(1,2,1)
hold on
% plot the bounds
plot([0,p.dt*nSteps],[0,0],'k:');
plot([0,p.dt*nSteps],[p.a,p.a],'r:');
plot([0,p.dt*nSteps],-[p.b,p.b],'r:');

stairs(t,y);
xlabel('RT (s)');
title('Unequal step size, equal probability')

%% 
% A second way to implement the diffusion model is to fix the up and down
% step sizes, and flip a biased coin to determine whether the walk goes up
% or down.

% The up and down step size is s*sqrt(dt) - same as the standard deviation
% of the step size above.
yStep = sqrt(p.dt)*p.s;
% The coin is biased by the drift rate.  Let's derive it.  If p is the
% probability of going up, then the expected dy of a given step is yStep*p
% - yStep*(1-p) = yStep*(2p-1) = s*sqrt(dt)*(2p-1). This should be equal to
% the drift rate, u*dt.  So s*sqrt(dt)*(2p-1) = u*dt.  Solving for p (prob)
% gives:
prob = .5*(sqrt(p.dt)*p.u/p.s + 1);

% To implement this, use the 'rand' and 'sign' functions:
dy = yStep*sign(rand(nSteps,nReps)-(1-prob));
y = cumsum(dy);
subplot(1,2,2)
hold on
% plot the bounds
plot([0,p.dt*nSteps],[0,0],'k:');
plot([0,p.dt*nSteps],[p.a,p.a],'r:');
plot([0,p.dt*nSteps],-[p.b,p.b],'r:');

stairs(t,y);
xlabel('RT (s)');
title('Equal step size, unequal probability')


%% Simulating RT distributions
%
% The walk should stop when the acculating variable hits one of the
% decision bounds.  The time step when this happens is the RT for that
% trial.

% We'll implement this with a loop over time.  We'll only keep track of the
% current values of y, and only increment the walks that haven't lead to a
% decision. 

nReps = 5000;  %number of simulated walks
p.dt = 1/1000; %time step (seconds)

%inialize some parameters:
alive = true(1,nReps);  %index vector for walks that haven't terminated
y = zeros(1,nReps); %Starting, and to be current position for each walk
response = zeros(1,nReps);  %will be filled with -1 or 1 for incorrect and correct
RT = zeros(1,nReps);  %will be filled with RT values in seconds

%for the coin-flip random walk, as before:
yStep = sqrt(p.dt)*p.s;
prob = .5*(sqrt(p.dt)*p.u/p.s + 1);

tStep = 0;
while sum(alive) %loop until all walks are terminated
    tStep = tStep+1;
    
    %for the 'randn' implementation, use this line: dy = p.u*p.dt +
    %p.s*sqrt(p.dt)*randn(1,sum(alive));

    %for the 'coin flip' implementation, use this line:
    dy =  yStep*sign(prob-rand(1,sum(alive)));

    %increment the 'living' walks
    y(alive) = y(alive)+ dy;

    %find the walks that reached the 'a' boundary
    aboveA = find(y>=p.a & alive);
    if ~isempty(aboveA)
        response(aboveA) = 1; %correct response
        alive(aboveA) = false; %'kill' the walk
        RT(aboveA) = tStep*p.dt; %record the RT
    end
    %find the walks that reached the '-b' boundary
    belowB = find(y<=-p.b & alive);
    if ~isempty(belowB)
        response(belowB) = -1;  %incorrect response
        alive(belowB) = false; %'kill' the walk
        RT(belowB) = tStep*p.dt; %record the RT
    end
end


%Plot histograms of simulated RTs

bins = linspace(0,max(RT),101); %Set time bins for histograms

figure(2)
clf
subplot(1,3,1) %correct
hist(RT(response==1),bins);
set(gca,'XLim',[min(bins),max(bins)]);
xlabel('RT (s)');
title('Correct');

subplot(1,3,2) %incorrect
hist(RT(response==-1),bins);
set(gca,'XLim',[min(bins),max(bins)]);
xlabel('RT (s)');
title('Incorrect');

subplot(1,3,3) %all
hist(RT,bins);
set(gca,'XLim',[min(bins),max(bins)]);
xlabel('RT (s)');
title('All');


%% Compare the simulated accuracy to the analytical solution. 

% Call the function 'expectedPC' which gives the expected probability
% correct.  

Psim = sum(response==1)/nReps;
Pexp = expectedPC(p);

disp(sprintf('simulated P(C): %%%5.2f, expected P(C): %%%5.2f',Pexp*100,Psim*100));


%% Compare the distribution of the simulated RTs to the analytical solution.

% Call the function 'expectedRTpdf' which gives the expected pdf 
% (histograms) of the RTs for correct answers.

RTpdfexp = (bins(2)-bins(1))*expectedRTpdf(p,bins,20);  %scale by time step (dt)

%Re-plot the histogram of correct responses
figure(3)
clf

n=hist(RT(response==1),bins);
n = n/sum(response==1);
stairs(bins,n,'b-');

% Draw the analytical solution on top.
hold on
plot(bins,RTpdfexp,'r-');
legend({'Simulated','Expected'});
xlabel('RT (s)');

%% to do:

% (1) Play around with the model parameters: drift rate, decision bounds - and
% see how the affect the accuracy and distribution of RTs.  Think about how
% a given data set (either behavioral or physiological) can help constrain
% these variables.  What is the signature of an increase in the decision
% bound p.b?  What is the signature for an increase in noise?  

% (2) Generate a psychometric function by calculating the percent correct
% as a function of drift rate, either through simulation or through the
% analytical solution.  See how the slope and threshold for this function
% varies with the mean and standard deviations of the drift rates.  Does
% this make sense? 




