% Chi2TestTutorial
%
% A simple tutorial of the chi-square distribution and chi-square test for whether
% a sample was drawn from a known distribution.
%
% Reference:  Statistics (Fifth edition)  William L. Hays, Harcourt Brace
%
% Requires: Statistics Toolbox
%
% 06/02/06      bx   Wrote it.
% 06/14/06      dhb  Editing and polishing.

%% Initialize
clear all; close all; 

%% What is the chi-square distribution?
%
% The chi-square distribution with N degrees of freedom is the sum of the
% squares of N draws from a N(0,1) distribution.
%
% Let's generate some draws from a chi-square and compare with the
% theoretical PDF.

% Calculate chi-squares, 1 df and compare with pdf returned by MATLAB
% function.  Need to multiply pdf by number of samples and histogram bin
% width to get the scale right.
nSamples = 5000;
r1 = randn(nSamples,1);
chisquare1 = r1.^2;
[hist1,x1] = hist(chisquare1,100);
y1 = nSamples*(x1(2)-x1(1))*chi2pdf(x1,1);
figure; clf;
subplot(1,2,1); hold on
bar(x1,hist1);
plot(x1,y1,'r','LineWidth',2);
title('Chi-square, 1 degree of freedom');

% Chi-squares, N df
df = 10;
r2 = randn(nSamples,df);
chisquareN = sum(r2.^2,2);
[histN,xN] = hist(chisquareN,100);
yN = nSamples*(xN(2)-xN(1))*chi2pdf(xN,df);
subplot(1,2,2); hold on
bar(xN,histN);
plot(xN,yN,'r','LineWidth',2);
title(sprintf('Chi-square, %d degrees of freedom',df));

%% Statistical Test
% Suppose we have a set of data and want to test the null hypothesis that
% this data was drawn from a known distribution.  This can be done with a
% chi-squared test.
%  
%   Null Hypothesis: Samples X are drawn according to known distribution R
%     H = 1 do no reject the null hypothesis.
%     H = 0 reject the null hypothesis. 
%
% The test works as follows.  You bin the observed data into some set
% number of bins.  This gives an empirical histogram.  The known
% distribution provides a prediction for how many samples should fall into
% each bin.  Then you compute the statistic
%
%  chi2 = sum_over_bins( (observedN-expectedN)^2/expectedN )
%
% Under the null hypothesis, and in the limit as the number of observations
% is large, this is distributed as chi-squared, with nBins-1 degrees of
% freedom.
%
% The following simulation uses the chi-squared test to check the null
% hypothesis that samples X are drawn from a distribution R.  The known
% distribution is taken as an exponential with mean muR.  The data are
% drawn from an exponential with mean muX.
%   Set muX = muR for the null hypothesis to be true.
%   Set them different for it to be false.
% You can then simulate out the test for various choices on number of bins, sample size,
% and alpha level.  (The alpha level is the probability that you'll reject
% the null hypothesis when it is in fact true.)

alpha = 0.05;       % Alpha value; significance level 
nSamples = 1000;    % Sample size (size of X)
nSimulation = 5000; % Number of simulation iterations
nBins = 5;          % Number bins to use in computing chi-squared statistic for test
muR = 1;            % Mean of R;
muX = 1;            % Mean of X;

% Choose bins based on a sample from R.  This is a pretty rough and ready
% way to choose the bins.  What is important is that the expected number of
% samples that will end up in any given bin not get too small, or you'll be
% dividing by a very small number when you compute the chi-squared test
% statistic.  I don't really know how small is bad.
binsc = [-Inf linspace(muR/nBins,2*muR,nBins-1) Inf];

% Find theoretical number of samples that should end up in each bin.  Code
% is sort of clunky to handle the special cases at each end.  Perhaps this
% could be done in a slicker way -- I didn't worry about it much.
clear expectedN
for i = 2:length(binsc)-1
    if (i == 2)
        expectedN(i-1) = nSamples*expcdf(binsc(i),muR);
    else
        expectedN(i-1) = nSamples*expcdf(binsc(i),muR)-sum(expectedN(1:i-2));
    end
end
expectedN(length(binsc)-1) = nSamples-sum(expectedN(1:length(binsc)-2));
expectedN = expectedN';
if (any(find(expectedN < 0.05)))
    fprintf('Some bins have small number of expected values, not a good idea\n');
end

% Make a plot of observed and expected number, for one draw from the
% expected distribution.
temp = exprnd(muR,nSamples,1);
n = histc(temp,binsc); n = n(1:end-1);
figure; clf; hold on
bar([binsc(2:end-1) 2*binsc(end-1)-binsc(end-2)],[n expectedN]);
xlabel('Bin value'); ylabel('Count'); title('Empirical and expected draws from R');

% Simulate many draws of size nSamples from X according to R, and compute the relevant
% chi-square statistic.
for i = 1:nSimulation  
    
    % Simulate X under null hypothesis.
    X = exprnd(muX,nSamples,1);
    observedN = histc(X,binsc); observedN = observedN(1:end-1);
    
    % Compute chi squared statistic.
    df = length(observedN)-1;
    chi2(i) = sum( (expectedN-observedN).^2./expectedN );

    % Calculate the p values and compare with alpha to make binary decision
    pval(i) = 1 - chi2cdf(chi2(i),df);
    H(i) = (pval(i)>=alpha);
end

% Check rejection rate
fprintf('Theoretical mean: %g, data mean: %g\n',muR,muX);
fprintf('Null hypothesis rejected %g%% of trials, alpha = %g%%\n',100*length(find(H == 0))/length(H),100*alpha);

% Look at distribution of chi2 statistic.  It should match the red curve
% only when muX = muR.
figure; clf; hold on
[hist1,x1] = hist(chi2,50);
y1 = nSimulation*(x1(2)-x1(1))*chi2pdf(x1,df);
bar(x1,hist1);
plot(x1,y1,'r','LineWidth',2);
title(sprintf('Chi-square test dist''n, %d degrees of freedom',df));



