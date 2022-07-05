%% Psychophysics Tutorial II: Threshold Estimation
%
% Written by G.M. Boynton for Psych 448 at the University of Washington
% Adapted for CSHL in June 2010, and revised for 2012 and 2014.
%
% The goal of this tutorial is to define the detection threshold. We'll
% simulate psychophysical data based on signal detection theory.  Then
% we'll fit the simulated data with a parametric psychometric function
% using a maximum likelihood procedure. Then we'll then pick off the
% stimulus intensity level that predicts 80% correct performance.
%
% It is suggested you look at the Psychophysics I tutorial first, which
% covers signal detection theory, ROC curves and the 2AFC paradigm.
%
% requires the folder 'Psychophysics' to be added to the path
%addpath('Psychophysics');

%% The psychometric function
%
% Real psychophysicists don't report proportion correct as a dependent
% measure.  They often don't report d-prime either.  Instead, they report
% the 'threshold', which is the stimulus intensity required to produce a
% desired level of performance.  A big advantage to the threshold is that
% it is in the physical units of the stimulus, which can be more easily
% compared across paradigms and labs.  Another advantage is that it can
% produce a family of stimuli that all presumably produce the same internal
% response.
%
% There are various ways to find the intensity that yields a criterion
% performance level.  We'll implement the a 3-down 1-up staircase method.
% In this tutorial we'll simulate a subject's performance by assuming a
% specific psychometric function that predicts performance in a 2AFC task.
% This function will be a Weibull function.

%%  A SDT model to generate psychometric functions.
%
% If you've gone through the Psycho_Tutorial_I tutorial you've learned that
% for a 2AFC experiment under signal detection theory (SDT), d-prime can
% be calculated from percent correct by:

PC = .8;
dPrime = sqrt(2)*norminv(PC)

% The inverse of this, calculating percent correct from d-prime is:

PC = normcdf(dPrime/sqrt(2))

% This means that if the physical stimulus strength is equivalent to
% d-prime, then the expected psychometric function should be a cumulative
% normal:

intensityList = linspace(0,4,101);
dPrime = intensityList;
PC = normcdf(dPrime,0,sqrt(2));
figure(1)
clf
plot(intensityList,PC,'b-');
ylabel('Proportion correct');
xlabel('Stimulus intensity (d-prime)');
title('This is a cumulative normal');
%%
% But if d-prime represents the signal-to-noise ratio for the internal
% response, it does not typically increase linearly with physical stimulus
% strength.  In most cases, like for luminance or contrast, it seems that
% d-prime increases in a compressive fashion with stimulus strength.  A
% simple compressive function is a power function with an exponent less
% than 1.
%
% If d-prime increases as a power function then then percent correct will
% no longer be a cumulative normal:

intensityList = linspace(0,30,101);
dPrime = intensityList.^.5;  %power function of stimulus intensity
PC = normcdf(dPrime/sqrt(2));
figure(1)
clf
plot(intensityList,PC,'b-');
ylabel('Proportion correct');
xlabel('Stimulus intensity');
title('This is not a cumulative normal.');

%%
% Plotting on a log-intensity axis:

figure(1)
clf
plot(log(intensityList),PC,'b-');
ylabel('Proportion correct');
xlabel('Stimulus intensity');
set(gca,'XTick',log([.5,1,2,4,8,16]));
logx2raw

hold on
title('This looks like a Weibull function');

%% The Weibull function
%
% The Weibull cumulative distribution function is one of a variety of
% functions that predict a psychometric function.  It has the general form:
%
% y = 1-exp(-(x/lambda)^k)
%
% where x is the stimulus intensity and y is the proportion correct.  Lambda
% and k are free parameters.  A reparameterized version for our purposes
% is:
%
% y = 1-(1-g)exp(-(kx/t)^b)
%
% where
%
% k = -log((1-a)/(1-g))^(1/b)
%
% With this parameterization, g is the performance expected at chance (0.5
% in our case for 2AFC), t is the threshold, and a is the performance level
% that defines the threshold.  That is, when x=t then y=a.  The parameter b
% determines the slope of the function.
%
% Here's our version of our function:

help Weibull

%%
% Let's plot a family of curves for a fixed threshold (p.b) and varying
% slope (p.t).  We'll plot the curves on a log-x axis.

intensityList = [.0625,.25,1,4,16,32];
x = exp(linspace(log(min(intensityList)),log(max(intensityList)),101));
p.t = 1;
bList = [.5,1,2,4];
figure(1)
clf
subplot(1,2,1)
y = zeros(length(bList),length(x));
for i=1:length(bList)
    p.b = bList(i);
    y(i,:) = Weibull(p,x);
end
plot(log(x),y')
set(gca,'XTick',log([.125,.25,.5,1,2,4,8,16]));
logx2raw;
legend(num2str(bList'),'Location','NorthWest');
xlabel('Intensity');
ylabel('Proportion Correct');
title('Varying b with t=0.3');

% Notice how that the slopes vary, but all curves pass through the same
% point.  This is where x=p.t and y =a (about 0.8)

%%
% Next we'll plot a family of curves for a fixed slope(p.t) and varying
% threshold (p.t).

p.b = 2;
tList = intensityList;
subplot(1,2,2)
y = zeros(length(tList),length(x));
for i=1:length(tList)
    p.t = tList(i);
    y(i,:) = Weibull(p,x);
end
plot(log(x),y')

legend(num2str(tList'),'Location','NorthWest');
xlabel('Intensity');
ylabel('Proportion Correct');
set(gca,'XTick',log(tList));
logx2raw
title('Varying t with b=2');

% See how for each curve, when the intensity is equal to the threshold
% (p.t), then the proportion correct is ~80%. So the variable p.t determines
% where on the x-axis the curve reaches 80%.
%
% By varying the two parameters, p.b and p.t, we can generate a whole
% family of curves that are good at describing psychometric functions in a
% 2AFC experiment.


%% The 3-down 1-up staircase
%
% The method of constant stimuli is simple but has the disadvantage that
% many of the trials are either too hard or too easy. These trials
% contribute relatively little to our estimate of the subject's threshold.
%
% Instead, we'll adjust the stimulus on each trial based on previous
% responses using a 3-down 1-up 'staircase' procedure. Whenever the subject
% makes an incorrect response, the next trial has an easier (larger)
% intensity.  If the subject gets three trials in a row correct, the next
% trial gets harder.
%
% Staircase procedures deserves a tutorial on its own.  A 3-down 1-up is
% the simplest.  There are many others, with names like Quest and Pest.

% First, we need to define the exponent for the compressive nonlinearity
% that transforms stimulus strength to d-prime. 0.3 is a good value to
% simulate the dimension of stimulus contrast.

p.noiseMean = 0;
p.sd = 1;
p.p = .3;

nTrials = 100;
currentIntensityNum = length(intensityList);  %start with the easiest intensity.
correctInARow = 0;
clear results
for i=1:nTrials
    results.intensity(i) = intensityList(currentIntensityNum);
    
    
    p.signalMean = results.intensity(i).^p.p; %compressive nonlinearity
    
    noise = randn(1)*p.sd + p.noiseMean;
    signal = randn(1)*p.sd + p.signalMean;
    
    results.response(i) = signal>noise;
    
    if results.response(i) == 0; %incorrect
        currentIntensityNum = currentIntensityNum+1; %make it harder
        correctInARow = 0;
    else
        correctInARow = correctInARow+1;
        if correctInARow == 3
            currentIntensityNum = currentIntensityNum-1; %make it easier
            correctInARow = 0;
        end
    end
    currentIntensityNum = min(max(currentIntensityNum,1),length(intensityList));
end

% Let's plot the results from this single staircase.

xTrial = 1:nTrials;
figure(2)
clf
stairs(xTrial,log(results.intensity),'b-');
set(gca,'YTick',log(intensityList));
logy2raw;
hold on
% and add red and green data points for correct and incorrect responses
%
plot(xTrial(results.response==0),log(results.intensity(results.response==0)),'bo','MarkerFaceColor','r');
plot(xTrial(results.response==1),log(results.intensity(results.response==1)),'bo','MarkerFaceColor','g');
xlabel('Trial number');
ylabel('Stimulus Intensity');

% See how the trials started out easy and how the stimulus intensity jumps
% down after three correct (green) responses in a row and jumps up after
% every incorrect response (red).  In this design, the size of the jump
% starts out large and decreases over the first 10 trials or so.
%
%%
%
% Then we compile the results to generate a psychometric function:

nCorrect = zeros(1,length(intensityList));
nTrials = zeros(1,length(intensityList));

for i=1:length(intensityList)
    id = results.intensity == intensityList(i);
    nTrials(i) = sum(id);
    nCorrect(i) = sum(results.response(id));
end

PC = nCorrect./nTrials;
sd = PC.*(1-PC)./sqrt(nTrials);  %pq/sqrt(n)

% Plot PC as a function of intensityList:

figure(3)
clf
errorbar(log(intensityList),PC,sd,'b','MarkerFaceColor','b','LineStyle','none');
hold on
plot(log(intensityList),PC,'bo','MarkerFaceColor','w');
set(gca,'XTick',log(intensityList));
set(gca,'YLim',[0,1.05])
xlabel('Stimulus Intensity');
ylabel('Proportion correct');
logx2raw
set(gca,'YLim',[0,1.05]);
title('3-down 1-up');


%% Threshold
%
% The goal is to measure the minimum stimulus intensity that is 'reliably'
% detected. But what does this mean? We could mean a stimulus intensity
% that corresponds to above-chance performance.  More commonly, we choose
% fix performance level, say 80% correct and estimate the intensity that
% produces this level of performance. This is important.  Beware of people
% who get excited about the ability to reliably 'see' stimuli below
% 'threshold'!
%
% With a 3-down 1-up staircase in the limit, the intensity level should
% gravitate toward a value that has an equal probability of getting easier
% and getting harder. Thus, the probability of getting 3 in a row correct
% is equal to 1/2.  This probability is the cube root of 1/2:

(1/2)^(1/3)

% So a 3-down 1-up staircase method will naturally adjust the intensity
% level to show trials that lead to about 80% correct.  Cool?  I think so.

% How do we estimate the stimulus intensity that produces 80% correct
% performance? Looking at the graph, if we interpolate, we can guess that
% it's about 0.08.  But linear interpolation is a bad idea.  First, it
% throws out information from all data except for two intensity values, and
% second it requires the measured psychometric function data to be
% monotonic.
%
% Instead, we'll get to the heart of this tutorial and fit a smooth curve
% to the psychometric function data and use this best-fitting curve to pick
% off the threshold.

%% A first guess at a set of parameters
%
% Let's choose some initial parameters for a slope and threshold. We can
% choose the values that we used to generate our data:

p.t = .25;
p.b = 1;
pGuess = p;

y= Weibull(pGuess,x);

hold on

figure(3)
plot(log(x),y,'k-');
logx2raw

% Is this a good fit? To answer this we need a measure of how well the
% curve fits the data. A common measure for curve fitting is a
% sums-of-squared error (SSE), which is the sum of the squared deviations
% between the curves and the data points.  However, for proportion-correct
% data like this, SSE is not appropriate because deviations along the
% proportion-correct do not have equal weights.  A 10% deviation for
% performance around 50% is less meaningful than a 10% deviation around
% 90%.  You can see this by looking at the error bars.  They're always
% small near 100% - pq/sqrt(n)...
%

%% Likelihood
%
% For proportion-correct data (or any data generated through a binary
% process), the appropriate measure is 'likelihood'. Here's how it's
% calculated;
%
% For a given trial, i, at a given stimulus intensity, xi, the Weibull
% function predicts the probability that the subject will get the answer
% correct. So if the subject did get respond correctly, then the
% probability of that happening is:
%
% pi = W(xi)
%
% Where xi is the intensity of the stimulus at trial i, and W(xi) is the
% Weibull function evaluated at that stimulus intensity.  The probability
% of an incorrect answer is, of course:
%
% qi = 1-W(xi)
%
% Assuming that all trials in the staircase are independent, then for a
% given Weibull function, the probability of observing the entire sequence
% of subject responses is:
%
% prod(pi) = prod(W(xi))
%
% for correct trials, and
%
% prod(qi) = prod(1-W(xi))
%
% for incorrect trials.
%
% If we let ri = 1 for correct trials and ri=0 for incorrect trials, like
% the vector results.response, then a simple way to calculate the
% probability of observing the whole sequence in one line is:
%
% $$\prod W(x_i)^{r_i}\left(1-W(x_i)\right)^{\left(1-r_i\right)}$$
%
% Here's how to do it in Matlab with our data set. First we evaluate the
% Weibull function for the stimulus intensityList used in the staircase:

y = Weibull(pGuess,results.intensity);

% Then we calculate the likelihood of observing our data set given this
% particular choice of parameters for the Weibull:

likelihood = prod( y.^results.response .* (1-y).^(1-results.response))

%%
% We want to choose values of p.t and p.b to make this number as large as
% possible. This is a really small number - sort of surprising since we
% thought that our choice of parameters for the Weibull was reasonably
% good.  The reason for this small number is that the product of a bunch of
% numbers less than 1 gets really small. So small that it gets into the
% tolerance for matlab to represent small numbers. To avoid this, we can
% take the logarithm of the equation above - this expands the range of
% small numbers so that we're not going to run into machine tolerance
% problems:

% The Matlab calculation of this is:
logLikelihood = sum(results.response.*log(y) + (1-results.response).*log(1-y))

%%
% NaN?  Turns out that this calculation can fail if the values of y reach 1
% because the log of 1-y can't be computed.  The way around this is to pull
% the values of y away from zero and 1 like this:

y = y*.99+.005;

%%
% This is called a 'correction for guessing' and it deals with the fact
% that subjects will never truly perform at 100% for any stimulus intensity
% because human subjects will always make non-perceptual errors, like motor
% errors, or have lapses in attention.
%
% Here's a re-calculation of the log likelihood:

logLikelihood = sum(results.response.*log(y) + (1-results.response).*log(1-y))


%%
% This is a more reasonable magnitude to deal with. It's negative because
% the log of a number less than 1 is negative. Larger numbers (less
% negative) still correspond to 'good' fits.
%
% Let's calculate the log likelihood for a different set of Weibull
% parameters:
%
% I've written a function called 'fitPsychometricFunction.m' that makes
% this calculation.  It takes in as arguments the structure 'p' holding the
% function's parameters, the results structure, and a string containing the
% name of the function to use ('Weibull') in our case.  It's almost exactly
% the same as the calculation above except that it reverses the sign
% (multiplies the log likelihood by -1) of the log likelihood.  This
% provides a positive number where smaller numbers mean better fits.  This
% is compatible with the optimization search algorithms that we'll discuss
% in a bit.
%
% Here's how to use it.
%
fitPsychometricFunction(pGuess,results,'Weibull')

%%
% Note that the value is exactly the negative of the previous calculation.
%
% Let's try a different fit:

pBetter.t =1;
pBetter.b = 1;
fitPsychometricFunction(pBetter,results,'Weibull')
%%
% We can visualize this by plotting this new prediction in red on the
% old graph:

figure(3)
y = Weibull(pBetter,x);
hold on
plot(log(x),y,'r-');

%%  The log-likelihood surface
%
% How do we find the parameters that maximize the log likelihood?  Matlab
% has a function that does this in a sophisticated way.  But for fun, let's
% just look at what the log likelihood is for a range of parameter values.
% Since there are two parameters for the Weibull function, we can think of
% the log likelihood as being a 'surface' with p.b and p.t being the x and
% y axes.
%
% List of b and t values:
bList = linspace(0,2,31);
tList = exp(linspace(log(.01),log(16),31));

% Loop through the lists in a nested loop to evaluate the log likelihood
% for all possible pairs of values of b and t.

logLikelihoodSurface = zeros(length(bList),length(tList));

for i=1:length(bList)
    for j=1:length(tList)
        pSamp.b = bList(i);
        pSamp.t = tList(j);
        y = Weibull(p,results.intensity);
        logLikelihoodSurface(i,j) = fitPsychometricFunction(pSamp,results,'Weibull');
    end
end

%%
% We can visualize the surface by using a contour plot

figure(4)
clf
contour(log(tList),bList,logLikelihoodSurface,50)
xlabel('t');
ylabel('b');
set(gca,'XTick',log(intensityList));
logx2raw
colorbar

%%
% We can plot the most recent set of parameters as a symbol on this contour
% plot.

hold on
plot(log(pGuess.t),pGuess.b,'o','MarkerFaceColor','k');
plot(log(pBetter.t),pBetter.b,'o','MarkerFaceColor','r');

% This figure is like a topographical map. Each line represents the values
% of b and t that produce an equal value for log likelihood. The fourth
% argument into 'contour' is a list of these values.
%
% The figure shows a series of concentric shapes with the smallest circling
% around t=0.12 and b = 1.6.  This is the 'bottom' of our surface and is
% close to the best choice of Weibull parameters.

%% Matlab's 'fminsearch' routine and 'fit.m'
%
% Searching through the entire grid of possible parameters is clearly an
% inefficient strategy (especially if there are even more parameters to
% deal with).  Fortunately there is a whole science behind finding the best
% parameters to minimize a function, and Matlab has incorporated some of
% the best in their function called 'fminsearch'.
%
% I find the way fminsearch is called a bit inconvenient. 'fit' allows us
% to put all variables for the function into a structure (feeding my
% obsession with structures), and it makes it easy to keep some parameters
% fixed and others to vary. So I've written a function that calls
% fminsearch called 'fit.m'. Here's how to use it:

help fit

%%
% The first argument into 'fit' is the name of the function to be
% minimized.  In our case it's 'FitPsychometricFunction'.  This function
% must have the specific form of the output being the value to be minimized
% and the first input argument being a structure containing the parameters
% that can vary.  The other input arguments can have any form.
%
% The next argument input into 'fit' is the starting set of parameters in a
% structure.  The third argument is a cell array containing a list of
% fields in this structure that you want to let vary.
%
% The remaining inputs are the same inputs that you'd put into the function
% to be minimized.
%
% Here's how to use it for the Weibull function:

% Starting parameters
pInit = p;

[pBest,logLikelihoodBest] = fit('fitPsychometricFunction',pInit,{'b','t'},results,'Weibull')

%%
% Note that the parameters changed from the initial values, and the
% resulting log likelihood is lower than the two used earlier (using
% 'pGuess and pBetter'.  Let's plot the best-fitting parameters on the
% contour plot in green:

plot(log(pBest.t),pBest.b,'o','MarkerFaceColor','b');

%%
% See how this best-fitting pair of parameters falls in the middle of the
% circular-shaped contour.  This is the lowest point on the surface.
%
% We now have the best fitting parameters.  Let's draw the best predictions
% in green on top of the psychometric function data:

y = Weibull(pBest,x);
figure(3)
plot(log(x),y,'b-');

%%
% By design, the parameter 't' is the stimulus intensity that should
% predict 80% correct performance.  We'll put the threshold in the title of
% the figure

title(sprintf('stimulus intensity thresold: %5.4f',pBest.t));

%%
% We can demonstrate this by drawing a horizontal line at ~80% correct to
% the curve and down to the threshold value on the x-axis:

plot(log([min(x),pBest.t,pBest.t]),(1/2)^(1/3)*[1,1,0],'k-');

%% Fixed and free parameters.
%
% Suppose for some reason we want to find the threshold while fixing the
% slope to be 3.  With 'fit' it's easy - all we do is change the list of
% 'free' parameters (the third argument sent to 'fit'):

[pBest,logLikelihoodBest] = fit('fitPsychometricFunction',pInit,{'t'},results,'Weibull')

%%
% We can show this on the contour plot.  The threshold will be the
% best-fitting value along the horizontal line where b=1:

figure(4)
plot(log([min(tList),max(tList)]),pBest.b*[1,1],'c-')
plot(log(pBest.t),pBest.b,'ko','MarkerFaceColor','c');

y = Weibull(pBest,x);
figure(3)
plot(log(x),y,'c-');

%% Increment thresholds
%
% So far we've only talked about detecting a signal against noise.  It's
% easy to use SDT to predict the ability to discriminate two stimuli that
% differ in strength.  An example of this is the increment threshold, which
% is the smallest increment in response that leads to some criterion
% performance.
%
% We can calculate the increment threshold function from the nonlinear
% transducer function using our SDT model.  Recall that in a 2AFC task, the
% threshold is the stimulus intensity (or increment) that produces a
% criterion level of performance.  In our definition, the level of
% performance is (1/2)^(1/3) = .7937 proportion correct.  Ths corresponds
% to a d-prime value of:

p.p = .3;
dPrime = sqrt(2)*norminv((1/2)^(1/3))

% Suppose we have the following baseline intensity values:
baseIntensityList = [0:10];

% We want to estimate the increments in stimulus intensity that produce a
% d-prime increment in internal response.  We'll find that increment using
% a binary search algorithm in which we refine our estimates within
% narrowing brackets of high and low estimates.

lo = zeros(size(baseIntensityList));  %list of increments that are too low
hi = 40*ones(size(baseIntensityList)); %these increments are too large

baseResponse = baseIntensityList.^p.p;
nIter = 20;
for i=1:nIter
    mid = (hi+lo)/2;  %increments midway between high and low
    %this is our best estimate of the increment in response
    incrementResponse = (baseIntensityList+mid).^p.p-baseResponse;
    %find which ones are too high and move the brackets accordingly
    tooHigh = incrementResponse>dPrime;
    hi(tooHigh) = mid(tooHigh);
    lo(~tooHigh) = mid(~tooHigh);
end

predIncrementThreshold = (hi+lo)/2;

hold on
figure(5)
clf
subplot(1,2,1)
plot(baseIntensityList,baseIntensityList.^p.p,'k.-');
set(gca,'XTick',baseIntensityList);
xlabel('Baseline intensity');
ylabel('Internal response');
title('Nonlinear transducer function');

set(gca,'YLim',[0,5]);
subplot(1,2,2)
plot(baseIntensityList,predIncrementThreshold,'k.-');
set(gca,'XTick',baseIntensityList);
xlabel('Baseline intensity');
ylabel('Increment threshold');
set(gca,'YLim',[0,40]);
title('Increment response function');

% This is Weber's law!  Increment thresholds increase with baseline
% contrast because of the compressive transducer function (power function)
% between stimulus intensity and internal response.  As the baseline
% intensity increases, the slope of the transducer function gets shallower
% so you need a larger and larger stimulus increment to get the signal mean
% to increase enough to obtain the d-prime needed for the threshold.

% Mess with the exponent (p.p) to see the relationship between the two
% curves.  

%%

% Let's simulate the increment thresholds for the same range of baseline
% stimulus strengths with our power-function model.  This will involve
% repeated 3-down 1-up staircases for varying baseline intensities.  

intensityList = [.0156,.0625,.25,1,4,16,32,64,128];  %for the staircase

nTrials = 1000;  %large number of trials

p.sd = 1;

incrementThreshold = zeros(size(baseIntensityList));
p.shutup = 'please';  %any field here suppresses the output to the command window for 'fit'
for j=1:length(baseIntensityList);  %loop through the baseline intensities
    p.noiseMean = baseIntensityList(j).^p.p; %set the noise mean for the staircase
    
    %start the staircase for this baseline intensity
    currentIntensityNum = length(intensityList);  %start with the easiest intensity.
    correctInARow = 0;
    clear results
    for i=1:nTrials
        results.intensity(i) = intensityList(currentIntensityNum);
        
        
        p.signalMean = (baseIntensityList(j)+results.intensity(i)).^p.p; %compressive nonlinearity
        
        noise = randn(1)*p.sd + p.noiseMean;
        signal = randn(1)*p.sd + p.signalMean;
        
        results.response(i) = signal>noise;
        
        if results.response(i) == 0; %incorrect
            currentIntensityNum = currentIntensityNum+1; %make it harder
            correctInARow = 0;
        else
            correctInARow = correctInARow+1;
            if correctInARow == 3
                currentIntensityNum = currentIntensityNum-1; %make it easier
                correctInARow = 0;
            end
        end
        currentIntensityNum = min(max(currentIntensityNum,1),length(intensityList));
    end
    nCorrect = zeros(1,length(intensityList));
    n = zeros(1,length(intensityList));
    for i=1:length(intensityList)
        id = results.intensity == intensityList(i);
        n(i) = sum(id);
        nCorrect(i) = sum(results.response(id));
    end
    
    PC = nCorrect./n;
    
    [pBest,logLikelihoodBest] = fit('fitPsychometricFunction',p,{'b','t'},results,'Weibull');
    incrementThreshold(j) = pBest.t;
end

figure(5)
subplot(1,2,2)
hold on
plot(baseIntensityList,incrementThreshold,'bo','MarkerFaceColor','w');

% There should be a reasonable match between our simulated data and the
% expected increment response function.  The fit should improve with
% increasing number of trials per staircase.

% What good is this? The goal of the psychophysicist is to take behavioral
% results (psycho) to learn something about the transformation between the
% physical stimulus (physics) and the internal response.  This simulation
% shows how you can predict in increment response function for a given
% nonlinear transducer function.  But we can go the other way too.  For a
% given psychophysical increment response function we can find the best
% nonlinear transducer function that predicts our data.  

% It turns out that for contrast increments, the psychophysical contrast
% increment threshold function (sometimes called threshold vs. contrast or
% TvC function) is surprisingly consistent with contrast response functions
% measured in V1 with fMRI.  Moreover, the effects of surround stimuli and
% attention affect psychophysical and fMRI responses in  a consistent way.



%% Exercises
%
% 1) Find your own stimulus intensity threshold by fitting the Weibull
% function to your own psychophysical data generated in Lesson 4.

% 2) Run 50 trials with a stimulus intensity level set at your own
% threshold and see how close the proportion correct is to 80%

% 3) What happens when you use different initial parameters when calling
% 'fit'?  Try some wild and crazy starting values.  You might find some
% interesting 'fits'.  This illustrates the need to choose sensible initial
% values.  If your function surface is convoluted then a bad set of initial
% parameters may settle into a 'local minimum'.
%
% 4) How would you obtain an estimate of the reliability of the threshold
% measure?  After all, p.t is just a statistic based on the values in
% 'results'.
%
% 5) Run the increment threshold simulation for an exponent of 1.  What
% happens and why?  



