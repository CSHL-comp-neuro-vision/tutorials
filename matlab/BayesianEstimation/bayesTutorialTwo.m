% coinflipBayesTutorial
%
% This is a simple example of Bayes estimation, that provides an approach
% to estimating the probability p that a coin will come up heads, after
% observing a series of flips of the coin.  
%
% 6/20/06   dhb  Wrote it, starting with some code that Barry Wark and I
%                wrote for something else.

% Revision: RK 6/2010

% Initialize workspace
clear; close all;

% Suppose someone is tossing a coin, and this coin has a probability p of
% coming up heads on each flip.  We don't know p.  We get a series of
% observations of whether the coin comes up heads or tails, and we want to
% use Bayes rule to update our estimate pEst of p.

% Define true value of p
actualP = 0.3;

% First step.  Define a prior on p.  To make things simple, we'll describe
% our prior on a discrete set of possibilities that p might take on.  We'll
% just start with a prior that is flat over the range 0-1.  You can play
% with different choices by setting this to something else
nPs = 100;                              % Number of bins over which we discretize p.
possiblePValues = linspace(0,1,nPs);    % The bins.  We know p is between 0 and 1.
priorPProbs = 1/nPs*ones(size(possiblePValues));
originalPriorPProbs = priorPProbs;

% Plot the prior
theFigure = figure; clf;
subplot(1,2,1); hold on
plot(possiblePValues,priorPProbs,'r');
xlabel('Possible Values of p');
ylabel('Probability');
title('Prior over p');
axis([0 1 0 0.5]); axis('square');

% Compute the prior mean, which we might take as our initial estimate of p.
priorMean = sum(priorPProbs .* possiblePValues);

% Loop over observations.  Each time, we'll obsever a random coin toss
% (driven by the true probability p specified above), compute the posterior
% using Bayes rule, and take the mean of the posterior to get our current
% estimate pEst of p.  Then we'll update the prior by substituting in the
% posterior, and repeat.
%
% The program updates a plot after each iteration, and the graph shows how
% it evolves.
%
% You can insert a "pause" command at the end of the loop if
% you want to step through by hand.
%
% As it marches through, it slows down because redrawing the figure takes
% more time with more observations.  If you want to start simulating out
% more trials, pull out the figures from the loop.

nObservations = 100;                     % Number of observations to simulate
observations = nan(nObservations, 1);
pEst = nan(nObservations, 1);
for i = 1:nObservations
    % Flip a coin to simulate this observation
    if (rand < actualP)
        observations(i) = 1;            % Heads
    else
        observations(i) = 0;            % Tails
    end
    
    % Compute posterior from prior and likelihood.  We loop over all the
    % possible values p can take on, and for each multiply the prior times
    % the likelihood.  This gives us the unnormalized posterior.  We
    % normalize at the end of the loop by ensuring that the posterior sums
    % to one.
    posteriorPProbs = nan(1, length(possiblePValues));
    for j = 1:length(possiblePValues)
        % The datum is a binary variable.  For each possible value of p, the
        % likelihood of observing a 1 is p, and the likelihood of observing
        % a zero is 1-p.
        if observations(i) == 1
            likelihood = possiblePValues(j);
        else
            likelihood = 1-possiblePValues(j);
        end

        % The unnormalized posterior is the prior times the likelihood.
        posteriorPProbs(j) = priorPProbs(j)*likelihood;
    end
    
    % Normalize the posterior so it sums to one.
    posteriorPProbs = posteriorPProbs/sum(posteriorPProbs);
    
    % Compute the estimate after observing the flip as the mean of the posterior.
    % You could also use the max, or some other function of the posterior.
    % But we'll leave exploring this as an exercise.
    pEst(i) = sum(posteriorPProbs .* possiblePValues);
    
    % Plot the posterior over the prior.  Prior in red, posterior in blue
    figure(theFigure); 
    subplot(1,2,1); hold off
    plot(possiblePValues,originalPriorPProbs,'r'); hold on
    plot(possiblePValues,posteriorPProbs,'b');
    xlabel('Possible Values of p');
    ylabel('Probability');
    title(sprintf('Posterior after %d flips',i));
    axis([0 1 0 0.5]); axis('square');
    
    % Plot the true value of p (red line) and our running estimates (blue circles).
    % The prior mean is shown as a red circle.
    % The green asterisks indicate whether each observation was a heads or a tail.
    figure(theFigure); 
    subplot(1,2,2); hold on
    plot(1:i,observations(1:i),'g*');
    plot(0,priorMean,'ro');
    plot(1:i,pEst(1:i),'bo');
    plot(1:i,pEst(1:i),'b');
    plot(0:i,actualP*ones(1,i+1),'r');
    axis([0 i 0 1]);
    xlabel('Observation');
    ylabel('p');
    title('Running estimates of p');
    axis('square');
    drawnow;
    
    % Make prior for next observation our current posterior.  Can you
    % convince yourself that this is a sensible move?  [Suppose I kept
    % the prior constant, waited until I had made all the observations, and
    % then computed the posterior.  Would I get the same answer as if I did
    % it incrementally as here?  Hint: I think it's only true if each coin flip is
    % independent.  You can write code that does it the other way if you want,
    % and compare.]
    priorPProbs = posteriorPProbs;
    
    % Uncomment this to wait for keypress each time through.
    % pause;
    
end

