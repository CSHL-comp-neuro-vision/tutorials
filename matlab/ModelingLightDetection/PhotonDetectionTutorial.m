%%  
% PhotonDetectionTutorial

addpath('dependencies_PhotonDetectionTutorial');

%%
% Concepts:
% photon detection
% Poisson Distributions
% Signal detection theory
% Baysian Estimation 
% Two alternative forced choice task

%   Author: Greg Field 2004
%   edited: GDF 2006 (pared down: believe it or not!) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The goal of this tutorial is to develop your intuitions
%  regarding Bayesian approaches to decision making.  A cost function
%  has not been included, but an interesting extension would be to add one.
%  This intuition will be built though the specific example
%  of detecting dim flashes of light in the visual system.  

%  We begin by discussing some basic properties of light and by describing
%  the signal and noise properties of rod
%  photoreceptor responses.  Familiarity with some basic properties of
%  Poisson distributions is assumed.  If you are unfamiliar with the
%  Poisson distribution there is an appendix at the end of the tutorial
%  that you may wish to run now.

%  Photons emitted by most light sources are independent of eachother.
%  Thus the probability that n photons are 
%  emitted by the light source in a time window, \delta t, is given by the 
%  Poisson equation P(n|N) = exp(-N)(N^n)/(n!), where N is the mean number
%  of photons emitted by the source, and n is a particular number of photons.  

%  Let's explore this idea with a little thought experiment.
%  Imagine we have a constant light source behind a shutter.  We
%  place a photomultiplyer tube (PMT) on the other side of the 
%  shutter.  A PMT is a device for counting photons much like
%  how a Geiger counter detects and counts radioactive particles.  

%  Now we quickly open and close the shutter for a given amount of time.
%  We repeat this experiment over and over, each time noting
%  how many photons are counted by the PMT.  Let's say that the
%  mean number of photons across all these experiments is 100.  
%  Incidently this is about the number of photons delivered to the 
%  cornea with a just detectable flash.  The distribution of counts
%  will look something like this...

NumSamples = 1000;
PoissonMean = 100;
Samples = poissrnd(PoissonMean, 1, NumSamples);

Bins = 50:0.5:150;
[HistPoiss, Histx] = hist(Samples, Bins);
bar(Histx,HistPoiss)
xlabel('samples')
ylabel('number')
title('Distribution of counts given by our PMT')

%  Notice these samples take on integer values.  We can't measure 91.3
%  photons on a particular trial.  Photons are not divisable.  

%  Now lets say that we want to attenuate the light by a 
%  factor of 100 (An operation similar to that performed by the optics of
%  the eye).
%  We could accomplish this by placing a
%  2 log unit neutral densitity filter in between
%  our light source and the PMT.  Now the mean number of 
%  counts will equal 1.  

NumSamples = 1000;
PoissonMean = 1;
Samples = poissrnd(PoissonMean, 1, NumSamples);

Bins = 0:0.5:10;
[HistPoiss, Histx] = hist(Samples, Bins);
bar(Histx,HistPoiss)
xlabel('samples')
ylabel('number')
title('Distribution of counts given by our PMT')

%  Notice that on many trials no photons are detected by our PMT.  This is
%  despite the fact that the mean flash strength was non-zero.
%  Now suppose we want to know the probability of our PMT
%  reporting a count of zero, one, or two photons on a given trial.  
%  This is very easy to calculate in MATLAB.

Probability_of_Zero_Photons = poisspdf(0, PoissonMean)
Probability_of_One_Photons = poisspdf(1, PoissonMean)
Probability_of_Two_Photons = poisspdf(2, PoissonMean)

%  Now imagine we need to decide whether a the shutter opened based on the
%  output of our PMT.  The probability of being wrong is equal to the
%  probability of detecting zero photons (0.36 at a mean flash strength of
%  1 photon).  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  We have a bit of intuition now about how light 
%  interacts with a detector.  We will now look at 
%  light responses in a rod photoreceptor and use
%  a Bayesian approach to detecting photons.

load SinglesData
load NoiseData

%  If we deliver a flash where the rod on average absorbs one photon, the
%  rod will sometimes absorb zero, one, two, etc., photons, just like the
%  PMT described above.  
%  Here are several examples of the response generated by a rod 
%  photoreceptor when it absorbes a single photon. 

figure(1)
clf
timepts = 0.001:0.001:1.901;
for cnt = 1:5
	subplot(5,1,cnt)
	plot(timepts, SinglesData(cnt,:))
    axis([0 1.5 -1 2])
	if cnt == 5
		xlabel('seconds')
	end
	if cnt == 3
		ylabel('normalized amplitude')
	end
	if cnt == 1
		title('Example single photon responses in a primate rod')
	end
end

%  In each of these examples the time of the flash was zero,
%  and the mean response has been normalized to have unit amplitude.

%  We can also look at the mean single photon response.

figure(2)
clf
plot(timepts, mean(SinglesData))
title('mean single photon response')
xlabel('normalized amplitude')
ylabel('seconds')
hold on

%  Let us also look at responses in the rod when no photons are absorbed. 

figure(1)
clf
for cnt = 1:5
	subplot(5,1,cnt)
	plot(timepts, NoiseData(cnt,:))
    axis([0 1.5 -1 2])
	if cnt == 5
		xlabel('seconds')
	end
	if cnt == 3
		ylabel('normalized amplitude')
	end
	if cnt == 1
		title('Examples where the rod failed to absorb a photon')
	end
end

%  Again we can also look at the mean zero photon response

figure(2)
plot(timepts, mean(NoiseData), 'k')
axis([0 2, -0.5 1])
title('mean noise event')
xlabel('normalized amplitude')
ylabel('seconds')
hold off

%  This ensemble of responses where the rod fails to absorb a photon
%  provides a characterization of the noise generated by the photoreceptor.
%  This noise is often called continuous dark noise because it consists of
%  continuous fluctuations in the current generated by the rod in darkness.
%  On some trials, the continuous dark noise can generate large noise
%  fluctuations that can easily be mistaken for single photon responses.  

%  We illustrate this by comparing how similar each response 
%  (to zero and one photon) is with the mean single photon response.
%  In mathematical terms, this is done by computing the dot product of
%  every response with the mean single photon response.  

%  Let's look at the distribution of rod noise.

MeanSingle = mean(SinglesData) - mean(NoiseData);
Template = MeanSingle ./ dot(MeanSingle, MeanSingle);
NoiseProjs = NoiseData * Template';
SinglesProjs = SinglesData * Template';
TotalProjs = [SinglesProjs; NoiseProjs];

figure(1)
clf
NumBins = 100;
[junkhist, Bins] = hist(TotalProjs, NumBins);

NoiseHist = hist(NoiseProjs, Bins);
stairs(Bins, NoiseHist, 'k') 
ylabel('number')
xlabel('response amplitude')
title('distribution of rod responses')
legend('noise distribution')
hold on

%  As we might have expected, this
%  distribution is centered on zero, and looks fairly Gaussian, although
%  it certainly has some skew off to negative values.  For our current 
%  purposes we will approximate this distribution of noise with a Gaussian.

% Fit noise with a normal
coef = [0.5];
global NumResponses
NumResponses = length(NoiseProjs);
NoiseCoef = nlinfit(Bins', NoiseHist', 'ZeroMeanGaussian', coef);
NoiseFit = ZeroMeanGaussian(NoiseCoef, Bins);
NoiseSD = NoiseCoef;

plot(Bins, NoiseFit, 'g')
legend('noise distribution', 'noise fit')

%  Now we will look at the distribution of single photon responses
%  projected along the mean single photon response.

SinglesHist = hist(SinglesProjs, Bins);
stairs(Bins, SinglesHist, 'r')
legend('noise distribution', 'noise fit', 'singles distribution')

% We will also fit this distribution with a Gaussian.
NumResponses = length(SinglesProjs);
global NoiseSD
NoiseSD = abs(NoiseCoef);
coef = [0.2,1]';
SinglesCoef = nlinfit(Bins', SinglesHist', 'SinglesDist', coef)

SinglesFit = SinglesDist(SinglesCoef, Bins);
SinglesMean = SinglesCoef(2);
SinglesSD = SinglesCoef(1);
plot(Bins, SinglesFit, 'b')
legend('noise distribution', 'noise fit', 'singles distribution', 'singles fit')
hold off

%  At this point we have two Gaussian distributions.  One describes 
%  the distribution of noise events, and the other describes the distribution
%  of singles photon responses.  (You may be worried about what happened
%  to multiphoton responses -- we will be working with flash strengths
%  where the probability of absorbing more than one photon is insignificant).
%  With these distributions we are ready
%  to think about the problem of deciding whether a given response
%  comes from the noise or signal distribution by applying some 
%  Bayesian ideas

%  Given these distributions, if we observed a 
%  response with amplidute 1.2, we would classify that response as coming
%  from (see plot of distributions) the single photon response 
%  distribution. If we saw a response with amplitude 0.2, 
%  it is more likely to be from the noise distribution.  

%  To generalize these ideas a bit we would look for the crossing point
%  between the singles distribution and the noise distribution.  Responses
%  falling to the right of that point are more likely to be singles, and 
%  those to the left are more likely noise.  In other words, compute the
%  probability N(r|0,\sigma_D) and N(r|1,\sigma_A), where N(r|\mu, \sigma)
%  is just the normal distirbution with mean \mu and a standard deviation
%  \sigma.  Classify r according to which distribution it is more likely 
%  to be sampled.  

%  In Bayesian
%  terminology we are using our "evidence" to dictate whether 
%  a response comes from the singles or noise distribution.  So
%  far we haven't said anything about priors.  Let's ignore them
%  for a moment more and see what happens with and without priors.

%  Imagine you are a retinal ganglion cell pooling over 5000 rod
%  photoreceptors in the periphery of the retina.  It is your job to report
%  whether you have seen a flash of light.  Let's say you use the strategy
%  outlined in the above paragraph to decide whether a flash was delivered
%  in your receptive field.  Let's generate a psychometric funtion for our
%  imaginary ganglion cell in the context of a two alternative forced choice 
%  task.  The alternatives from which you choose are whether a flash
%  was or was not delivered.  The bit of code below runs a simulation that does
%  this.  Don't worry about the details of what it is doing, I will biefly
%  explain that at the end of the tutorial.

%  This bit of code will take a couple minutes to run -- get some coffee,
%  or better yet, a beer.
%  This bit intitializes some parameters
FlashStrengths = [0.0007 0.001 0.003 0.007 0.01 0.03 0.07];
NumFlashes = length(FlashStrengths);
RFSize = 5000;
NumTrials = 25;
Verbose = 0;
ThermalRate = 0;
NumPhotons = 10;
%  This bit simulates rod responses and calculates probability correct
load ThermalDist
for flash = 1:NumFlashes
    fprintf('processing flash %d photons per rod \n', FlashStrengths(flash))
    [SimNoise, SimSingles] = GenerateCombResponses(SinglesSD, NoiseSD, SinglesMean, 0, FlashStrengths(flash), RFSize, NumTrials, ThermalDist, ThermalRate, Verbose); 
    TempLin = LinearCombDiscrim(SimSingles, SimNoise, SinglesMean, SinglesSD, 0, NoiseSD, FlashStrengths(flash), NumPhotons, Verbose);
    LinPCorrect(flash) = TempLin;
end
%  Fit and plot the simulated ganglion cell responses.
coef = 0.003;
fitcoef = nlinfit(FlashStrengths', LinPCorrect', 'cumulative_gaussian', coef);
LinThreshold = norminv(0.75, 0, fitcoef);  
semilogx(FlashStrengths, LinPCorrect, 'ko')
hold on
semilogx(FlashStrengths, cumulative_gaussian(fitcoef, FlashStrengths), 'k')
xlabel('flash strength')
ylabel('probability correct')
hold off

%  This graph looks like your standard psychometric function.  We are 
%  plotting probability correct as a function of flash strength, and
%  we observe a standard sigmoidal relation. Let's compare this 
%  performance to that of an array of optimal photon
%  detectors (PMT's without any noise).  

OptFlashStrengths = [0.00003 0.00007 0.0001 0.0003 0.0007 0.001];
OptimalDetector = 1 - (exp(-OptFlashStrengths * RFSize) / 2);
hold on
semilogx(OptFlashStrengths, OptimalDetector, 'r');
hold off

%  This red line is the theoritically optimal performance.  
%  The threshold of the optimal photon detector is about 2 log units
%  lower than that of our pool of rod photoreceptors. Why is there such a
%  discrepancy between the performance of the optimal photon detector and
%  the pool of rod responses.  After all, when we plotted the distribution
%  of rod signals along side the distribution of noise, there was not a lot
%  of overlap between these distributions. 

%  The answer lies in the fact that we have not taken a complete Bayesian approach
%  to deciding whether photons have been absorbed in the receptive field
%  of our imaginary ganglion cell.  Essentially we have NOT weighted the
%  evidence by the priors.  If you know the flash strength, then you know
%  the relative probability of getting a response from the noise distribution
%  versus absorbing a photon (this just comes from Poisson statistics).
%  Thus, the priors are given by the flash strength.  
%  Let's replot the distribution of noise and single photon
%  responses and try to generate some intuition about how applying
%  priors to the evidence changes how we classify responses.

figure(1)
clf
subplot(2,1,1)
plot(Bins, 0.5 * normpdf(Bins, SinglesMean, sqrt(SinglesSD^2 + NoiseSD^2)), 'k', Bins, 0.5 * normpdf(Bins, 0, NoiseSD), 'r')
xlabel('response amplitude')
ylabel('probability density')
legend('singles dist', 'noise dist')

%  In this picture the probability of absorbing a photon versus not
%  absorbing a photon is roughly equal.  This is true for flash strengths
%  around 1 photoisomerization.  Imagine the flash strength is 0.01
%  photoisomerizations.  What happens to these two distributions?

subplot(2,1,2)
plot(Bins, 0.01 * normpdf(Bins, SinglesMean, sqrt(SinglesSD^2 + NoiseSD^2)), 'k', Bins, 0.99 * normpdf(Bins, 0, NoiseSD), 'r')
xlabel('response amplitude')
ylabel('probability density')
legend('singles dist', 'noise dist')

%  Now the singles distribution has so little probability mass that you can
%  barely see it.  Let's replot this with the abscissa in log units.

figure(2)
clf
subplot(2,1,1)
semilogy(Bins, 0.5 * normpdf(Bins, SinglesMean, sqrt(SinglesSD^2 + NoiseSD^2)), 'k', Bins, 0.5 * normpdf(Bins, 0, NoiseSD), 'r')
title('flash strength ~1 photons / rod')
xlabel('response amplitude')
ylabel('probability density')
subplot(2,1,2)
semilogy(Bins, 0.01 * normpdf(Bins, SinglesMean, sqrt(SinglesSD^2 + NoiseSD^2)), 'k', Bins, 0.99 * normpdf(Bins, 0, NoiseSD), 'r')
title('flash strength ~0.01 photons / rod')
xlabel('response amplitude')
ylabel('probability density')
legend('singles dist', 'noise dist')

%  The top panel shows the singles and noise distribution for a flash
%  strength of 1 photoisomerization. The bottom panel shows the
%  distributions for a flash strength of 0.01 photoisomerizations.  Notice
%  that the crossing point between these two distributions has moved
%  substantially to higher values in the bottom panel.  At even smaller
%  flash strengths, this crossing point will move to yet higher values.  

%  This implies that as we move to lower and lower flash strengths, our
%  criteria for deciding wether a response is generated by the absorption
%  of a single photon versus a noise event should increase.  
%  This is accomplished naturally in a Bayesian framework.
%  We imploy this framework in the following section.

%  Re-simulate flash responses...  This will take a couple minutes.
FlashStrengths = [0.00007 0.0001 0.0003 0.0007 0.001 0.003 0.007 0.01 0.03 0.07];
NumFlashes = length(FlashStrengths);
RFSize = 5000;
NumTrials = 25;
Verbose = 0;
ThermalRate = 0;
NumPhotons = 10;
for flash = 1:NumFlashes
    fprintf('processing flash %d photons per rod \n', FlashStrengths(flash))
    [SimNoise, SimSingles] = GenerateCombResponses(SinglesSD, NoiseSD, SinglesMean, 0, FlashStrengths(flash), RFSize, NumTrials, ThermalDist, ThermalRate, Verbose); 
    TempOpt = OptimalCombDiscrim(SimSingles, SimNoise, SinglesMean, SinglesSD, 0, NoiseSD, FlashStrengths(flash), NumPhotons, Verbose);
    TempLin = LinearCombDiscrim(SimSingles, SimNoise, SinglesMean, SinglesSD, 0, NoiseSD, FlashStrengths(flash), NumPhotons, Verbose);

    OptPCorrect(flash) = TempOpt;
    LinPCorrect(flash) = TempLin;
end
coef = 0.003;
Linfitcoef = nlinfit(FlashStrengths', LinPCorrect', 'cumulative_gaussian', coef);
LinThreshold = norminv(0.75, 0, Linfitcoef); 
coef = 0.0007;
Optfitcoef = nlinfit(FlashStrengths', OptPCorrect', 'cumulative_gaussian', coef);
OptThreshold = norminv(0.75, 0, Optfitcoef);

% plot performance and fits from simulation
clf
semilogx(FlashStrengths, LinPCorrect, 'ko')
hold on
semilogy(FlashStrengths, OptPCorrect, 'bo')
semilogx(FlashStrengths, cumulative_gaussian(Linfitcoef, FlashStrengths), 'k')
semilogx(FlashStrengths, cumulative_gaussian(Optfitcoef, FlashStrengths), 'k')
xlabel('flash strength')
ylabel('probability correct')
hold off

%  Let's compare this performance to that of an array of optimal photon
%  detectors (PMT's without any noise).  
OptFlashStrengths = [0.00003 0.00007 0.0001 0.0003 0.0007 0.001];
OptimalDetector = 1 - (exp(-OptFlashStrengths * RFSize) / 2);
hold on
semilogx(OptFlashStrengths, OptimalDetector, 'r');
hold off

%  There is a very large improvement in performance when the priors
%  are used to weight the evidence.  This shows the importance
%  of considering the priors in a detection task.  Quite remarkably
%  the retina appears to be optimized for detecting light
%  at absolute threshold, and weights rod responses according
%  to a prior where only 1 in 10,000 rods is absorbing a photon
%  (Field and Rieke 2002).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Brief Intro to Poisson Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Here is a very brief introduction to Poisson statistics.  
%  If you have not already encountered this distribution, 
%  you will use it often in the context of this course and beyond.  
%  In a nutshell, the Poisson distribution describes the 
%  statistics of independent events.  Perhaps the most important 
%  property is that the variance of a Poisson distribution 
%  equals the mean.  The mathematical form of the distribution
%  is P(x) = exp(-n)(n^x)/(x!).  Let's look at a few different 
%  Poisson distributions.

NumSamples = 100;
PoissonMean = 1;
Samples = poissrnd(PoissonMean, 1, NumSamples);

NumBins = 20;
Bins = 0:0.5:(0.5 * NumBins);
[HistPoiss, Histx] = hist(Samples, Bins);
bar(Histx,HistPoiss)
xlabel('samples')
ylabel('number')
title('Poisson distribution')

%  The first thing to notice about this distribution 
%  is that samples take on only discrete values, e.g.,
%  0, 1, 2,...  Notice also that we set PoissonMean = 1,
%  and the mean of our samples is very nearly one

mean(Samples)

%  Also recall that the variance of a Poisson 
%  distribution equals the mean.

var(Samples)

%  What happens when the mean of the Poisson distribution
% is changed?  Here are some examples.

PoissonMean = [0.1 0.5 1, 2, 4];
NumMeans = length(PoissonMean);
for cnt = 1:NumMeans
	TempSamples = poissrnd(PoissonMean(cnt), 1, NumSamples);
	[HistPoiss, Histx] = hist(TempSamples, Bins);
	subplot(NumMeans, 1, cnt)
	bar(Histx,HistPoiss)
	xlabel('samples')
	ylabel('number')
end

%  In each distribution shown above the samples are again integers.
%  This is true even when the mean is not an integer, such as the
%  top two plots.  You can also check that the variance equals
%  the mean in each of the distributions above.