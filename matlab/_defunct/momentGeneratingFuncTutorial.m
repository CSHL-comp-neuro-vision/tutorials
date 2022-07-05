% momentGeneratingFuncTutorial
% Michael Shadlen 2003
% Revision date June 28, 2004

addpath('momentGeneratingFunc');

% This tutorial introduces the moment generating function. It complements (and
% repeats to some extent) the section in my mathematica tutorial on Wald's
% identity (Wald_identity.nb). 
%
% My goal is to develop a few basic intuitions that are essential for
% understanding Wald's identity and its connection to the psychometric function.
% Here are the topics you need to understand.
% The basic definition of the MGF
% The use of the MGF to calculate moments
% The MFG for a random variable, z, that is the sum of two independent random
% variables, x + y, is the product of the MGF's associated with these two random
% variables. 
% Intuitions for why the moment generating function looks the way it does. One
% of the more important goals is to gain an intuition for how changes in the
% properties of a random variable, its mean and variance, affect a special point
% on the moment generating function where it reaches 1: the special root,
% theta1. The reason we care about this is spelled out in the mathematica
% tutorial. It is the connection between Wald's Identity and the psychometric
% function.
%

% This short tutorial accomplishes two things. In part 1, we define the moment
% generating function (mgf) and use it to compute moments (hence the name). In
% the second part, we look for the existence of a special root to the mgf. This
% is the nonzero solution of the equation mgf(z) = 0. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1. The Moment Generating Function. Definition and properties
% 
% The MGF is a transform of a probability function or a probability density
% function. The PDF, f(x),  is the likelihood (or probability) of observing a
% random value, x. The MFG, M(theta) is a function of a new variable. Think of
% it the same way you would think of a fourier or laplace transform.

% Let's begin with an example of something familiar.
clear all

% define an x axis. We use x for the random variable
dx = .01;            % use this to control the sampling interval for x
x = [-10:dx:10]';

% Let's start with the normal distribution. We'll play with other distributions
% later.
u = 1; s = 1.2;       % mean and standard deviation
pdf = normpdf(x,u,s);

% Look at the probability density function (pdf)
figure(1), clf, hold on, plot(x,pdf)
axprefs(gca)
xlabel('x'); ylabel('Probability density')

% The moment generating function is a function of a new variable, theta. It is
% the expected value of exp(theta*x), where x is our random variable. Notice
% that for any value of theta, the expectation of e^(theta*x) is a number.

% define a theta axis. You may have to try a different range, depending on your
% choice of mean, stdev, etc.
dtheta = .01;       %  the grain of the axis
theta = [-3:dtheta:1]';

% The mgf at each value of theta is the expectation of exp(thata(i) * x)
% The expectation is just the  sum of exp(thata(i) * x) with each term weighted
% by the probability of of observing x. For example, when theta = 0
i = find(theta==0)  % i is the index to the theta vector, such that theta(i) is 0.
y = trapz(x, exp(theta(i).*x)  .* pdf)
% Convince yourself that all moment generating functions must return 1 at M(0).

% Try some other values.

% I wrote a function to handle processing of a vector of values in theta
y = mgf(theta,x,pdf);
figure(2), clf, hold on, plot(theta,y)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');

% The goal of the next section is to get a feel for why the MGF looks the way it
% does. Let's just notice a few things about this function. We already noticed
% that the value at theta=0 is 1. It should also be obvious that all values must
% be greater than or equal to 0.  Remember, it's just a weighted average of
% exp(...).  Notice that on one side of 0, the function keeps getting bigger. On
% the other side of 0, it dips below 1 and then rises again, crossing 1 at
% another point. That crossing point  turns out to be very important for
% psychophysics. We'll spend time on it later. 
 % For the moment (no pun intended) let's develop an
% intuition for these properties.

% The moment generating function has its name because if you take its derivative
% at theta=0, you get the 1st moment of x. In other words the mean. If you take
% the second derivative, again evaluated at theta=0, you obtain the 2nd moment:
% the average of the square of x. And so forth for higher derivatives and
% moments. This is easy to see in the math. Look at section 1.2 of
% Wald_identity.nb.  Here, let's convince ourselves of this using numerical
% approximations.

% The first moment is the expectation of x (i.e., the mean). We know what to
% expect, because we set the mean = u
u
% Or we can get this from 
normstat(u,s)
% or by calculating the expectation as the weighted sum of x's (weights from the
% pdf)
trapz(x,x .* pdf)
% OR we can differentiate the mgf and evaluate the derivative at 0
% To differentiate, we use diff
d = diff(y) ./ dtheta;
% These are right sided differences x(2)-x(1), x(3)-x(2), ... 
% the derivative at theta=0 is not represented directly. We'll approximate it by
% looking at the two vals on either side of 0
n = find(theta==0)
% interpolate to approximate the value of the derivative at 0.
derivAt0 = mean([d(n-1) d(n)])

% This also provides some intution for why the MGF's appearance. At theta=0, we
% know that the tangent to the curve should have slope = u.

% show the deriv on the fig
L = abs(theta) < .35;
tanLine = derivAt0*theta +1;
plot(theta(L), tanLine(L),'k')

% The second moment is the expectation of x^2. We also know what this should be,
% because we set u and s, above.
% Since the variance is <x^2> - <x>^2,  the second moment should be
m2 = s^2 + u^2
% or calculating directly
trapz(x,x.^2 .* pdf)
% or using the moments
dd = diff(y,2) ./ dtheta^2;
dd(n-1)         % 2nd derivative evaluated at 0

% Looking at the curve, we're not surprised that it is convex up at theta=0. The
% second moment has to be positive. If the variance were larger, then the degree
% of convexity would increase, making for a tighter U-shaped curve. Let's pay
% attention to the point on the left portion of the curve where it crosses the
% dashed line. The value of theta where this occurs is theta_1. This value is
% going to turn out to be very important to us. So pay close attention to this
% point in the following graphs.
%
u1 = 1; s1 = 1.2;       % mean and standard deviation
u2 = u1; s2 = 1.5 * s1;  % same mean, bigger stdev
pdf1 = normpdf(x,u1,s1);
pdf2 = normpdf(x,u2,s2);
MGF1 = mgf(theta,x,pdf1);
MGF2 = mgf(theta,x,pdf2);
% Look at the probability density function (pdf)
figure(1), clf, hold on, plot(x,pdf1, x, pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')
title('same mean, larger stdev')
% and MGF
figure(2), clf, hold on, plot(theta,MGF1,theta, MGF2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');
L = abs(theta) < .35;
tanLine1 = u1*theta +1;
tanLine2 = u2*theta +1;
plot(theta(L), tanLine1(L),theta(L), tanLine2(L))
title('same mean, larger stdev (green)')

% This is an important intuition to hang on to. If the convexity were greater,
% the second point where the curve crosses the horizontal line would move closer
% to theta=0.  In other words theta_1 would be smaller in absolute magnitude. 

% If the variance were smaller, the convexity would be lower and theta_1 would
% move off to the left, further from 0.
%
u1 = 1; s1 = 1.2;       % mean and standard deviation
u2 = u1; s2 = .7 * s1;  % same mean, smaller stdev
pdf1 = normpdf(x,u1,s1);
pdf2 = normpdf(x,u2,s2);
MGF1 = mgf(theta,x,pdf1);
MGF2 = mgf(theta,x,pdf2);
% Look at the probability density function (pdf)
figure(3), clf, hold on, plot(x,pdf1, x, pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')
title('same mean, smaller stdev')
% and MGF
figure(4), clf, hold on, plot(theta,MGF1,theta, MGF2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');
L = abs(theta) < .35;
tanLine1 = u1*theta +1;
tanLine2 = u2*theta +1;
plot(theta(L), tanLine1(L),theta(L), tanLine2(L))
title('same mean, smaller stdev')



% If the mean were larger, but the convexity did not change,
% then the slope at 0 would increase, and that would push theta_1 away from 0.
%
u1 = 1; s1 = 1.2;       % mean and standard deviation
u2 = 1.5*u1; s2 = s1;  % larger mean, same stdev
pdf1 = normpdf(x,u1,s1);
pdf2 = normpdf(x,u2,s2);
MGF1 = mgf(theta,x,pdf1);
MGF2 = mgf(theta,x,pdf2);
% Look at the probability density function (pdf)
figure(5), clf, hold on, plot(x,pdf1, x, pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')
title('larger mean, same stdev')
% and MGF
figure(6), clf, hold on, plot(theta,MGF1,theta, MGF2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');
L = abs(theta) < .35;
tanLine1 = u1*theta +1;
tanLine2 = u2*theta +1;
plot(theta(L), tanLine1(L),theta(L), tanLine2(L))
title('larger mean, same stdev')

% If the mean were smaller, but the convexity did not change,
% then the slope at 0 would decrease, and that pull theta_1 closer to the origin
%
u1 = 1; s1 = 1.2;       % mean and standard deviation
u2 = .75*u1; s2 = s1;  % smaller mean, same stdev
pdf1 = normpdf(x,u1,s1);
pdf2 = normpdf(x,u2,s2);
MGF1 = mgf(theta,x,pdf1);
MGF2 = mgf(theta,x,pdf2);
% Look at the probability density function (pdf)
figure(7), clf, hold on, plot(x,pdf1, x, pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')
title('smaller mean, same stdev')
% and MGF
figure(8), clf, hold on, plot(theta,MGF1,theta, MGF2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');
L = abs(theta) < .35;
tanLine1 = u1*theta +1;
tanLine2 = u2*theta +1;
plot(theta(L), tanLine1(L),theta(L), tanLine2(L))
title('smaller mean, same stdev')

% Now check this out. We can pit these two tendencies against each other. It
% turns out that their effects cancel if we change the mean and variance by the
% same amount. 

u1 = 1; s1 = 1.2;       % mean and standard deviation
u2 = 1.5*u1; s2 = sqrt(1.5)*s1;  % proportionately larger mean and variance
pdf1 = normpdf(x,u1,s1);
pdf2 = normpdf(x,u2,s2);
MGF1 = mgf(theta,x,pdf1);
MGF2 = mgf(theta,x,pdf2);
% Look at the probability density function (pdf)
figure(9), clf, hold on, plot(x,pdf1, x, pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')
title('mean and variance scaled identically')
% and MGF
figure(10), clf, hold on, plot(theta,MGF1,theta, MGF2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');
L = abs(theta) < .35;
tanLine1 = u1*theta +1;
tanLine2 = u2*theta +1;
plot(theta(L), tanLine1(L),theta(L), tanLine2(L))
title('mean and variance scaled identically')

% That's an important observation. For the gaussian distribution this special crossing
% point, theta_1, is a function of the ratio of variance and mean. If they are
% kept the same, theta_1 does not change. That's not to say that the MGF doesn't
% change. You can see from Figure 10 that it does. But that value, theta_1, is
% going to turn out to be pivotal to the performance.

% Already, you should be scratching your head skeptically. Shouldn't performance
% have something to do with signal to noise ratio, hence mean and standard
% deviation? Yet the theta_1 depends on mean and *variance*.  Store the thought.
% We have yet to connect theta_1 to the psychometric function.


% Let's expand our intuitions by looking at MGFs associated with mean=0 or mean
% < 0.  Let's start with the former.

% When the mean is zero, we know that the slope of the MGF at theta=0 is flat.
% Moreover, the function is convex up. So there's no 2nd crossing of the
% horizontal line. theta_1 is effectively 0. 
%
u1 = 0; s1 = 1.2;       % mean 0; same standard deviation as above
u2 = 0; s2 = sqrt(1.5)*s1;  % same mean, bigger variance
pdf1 = normpdf(x,u1,s1);
pdf2 = normpdf(x,u2,s2);
MGF1 = mgf(theta,x,pdf1);
MGF2 = mgf(theta,x,pdf2);
% Look at the probability density function (pdf)
figure(11), clf, hold on, plot(x,pdf1, x, pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')
title('mean = 0, green has larger variance')
% and MGF
figure(12), clf, hold on, plot(theta,MGF1,theta, MGF2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');
L = abs(theta) < .35;
tanLine1 = u1*theta +1;
tanLine2 = u2*theta +1;
plot(theta(L), tanLine1(L),theta(L), tanLine2(L))
title('mean = 0, green has larger variance')

% When the mean is less than 0, everything simply flips around the y-axis
% I think this only seems obvious, but if you think about it (or do the math)
% you'll see that  if pdf2(x) = pdf1(-x), then the weighted averages of
% exp(theta*x) with weights given by pdf2 have to correspond to the opposite
% signed theta. In any case we can see this graphically.


% I'll extend the range of theta so we can appreciate the symmetry.
th2 = [-3:dtheta:3]';

u1 = 1; s1 = 1.2;       % mean and standard deviation
u2 = -u1; s2 = s1;  % opposite sign mean, same stdev
pdf1 = normpdf(x,u1,s1);
pdf2 = normpdf(x,u2,s2);
MGF1 = mgf(th2,x,pdf1);
MGF2 = mgf(th2,x,pdf2);
% Look at the probability density function (pdf)
figure(13), clf, hold on, plot(x,pdf1, x, pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')
title('opposite means')
% and MGF
figure(14), clf, hold on, plot(th2,MGF1,th2, MGF2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([-.1 min(get(gca,'YLim'))]) 5])
line(get(gca,'XLim'),[1 1], 'LineStyle','--');
L = abs(th2) < .35;
tanLine1 = u1*th2 +1;
tanLine2 = u2*th2 +1;
plot(th2(L), tanLine1(L),th2(L), tanLine2(L))
title('opposite means')

% The ability to think about the negative version of a random variable will be
% important in a moment when we consider sums of RVs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The MGF of a sum of RVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% I may fill this in with more intuition later. But here's the bottom line.
% Suppose x1 has distribution pdf1(x) and moment generating function M1(theta),
% and suppose x2 has distribution pdf2(x) and moment generating function
% M2(theta).  Consider the random variable, z. It has a pdf that is the
% convolution of pdf1(x) and pdf2(x). It has a moment generating function that
% is the product M1(theta).*M2(theta). 
%
% It may seem surprising that the pdf(z) is a convolution, but think about it.
% Suppose we know the value, x1. Then the conditional probability of observing
% any z as the sum x1 + x2  is the probability of choosing x2 = z-x1. That's pdf2(z-x1). To get the
% probability of getting z from the any sum, it's just a matter of integrating
% this conditional probability across all possible values of x1. 
%  
% pdf(z) =  integral of pdf1(x1) .* pdf2(z-x1)
%
% In other words pdf(z) = conv(pdf1, pdf2)
% 
% I'm not writing real matlab code here. We would have to be careful about axes.
%
% I'm not going to go through the math, but instead I'm appealing to an intution
% that I hope we share. If we were dealing with functions of time or space, we
% would know that the fourier transform of pdf(z) would be the product of the
% fourier transforms of pdf1 and pdf2. The moment generating function is a lot
% like a fourier transform.
%
% I may flesh this out one day.
%
% There are two points relevant to our topic. First, if we were to make a new RV
% from the difference of two RVs, we need simply multiply the M1(theta) by
% M2(-theta). That's because we're adding x1 and -x2. Second, when we look at
% Wald's identity, we should not be surprised the the moment generating function
% for a sum of RVs is equal to the product of MGFs for the increments. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2. The mfg of a variety of RVs formed by taking the difference of two RVs
% has a nonzero theta_1. Recall that if z is the sum of two RVs, x + y, each
% with mfg Mx(theta) and My(theta), the mgf associated with z is
%
% Mz = Mx(theta).*My(theta). 
% Subtraction is like adding one of the RVs with its sign switched. This is
% Mx(theta).*My(-theta)
% 

% Let's look at the dfference between two RVs described by two weibull
% distributions. 
%
x = [0:dx:20];
%
pdf1 = wblpdf(x,.02,3);
pdf2 = wblpdf(x,.01,3);
%
figure(1), clf
plot(x,pdf1,x,pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')

theta = [-1:dtheta:1]';
% MGF for the difference of the 2nd Weibull minus the 1st. Notice the
% negative theta on the 1st mgf.
y = mgf(-theta,x,pdf1) .* mgf(theta,x,pdf2) ;

figure(2), clf, hold on
plot(theta,y,'k-','LineWidth',2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([.5 .5*min(y)]) 2])
hl = line(get(gca,'XLim'),[1 1])
set(hl,'LineStyle','--','Color','k')

% Use the MGF to get the mean difference. 
d = diff(y) ./ diff(theta);
n = find(theta==0);
% interpolate to approximate the value of the derivative at 0.
derivAt0 = mean([d(n-1) d(n)]);
meanEst = derivAt0
% show the deriv on the fig
L = abs(theta) < .2;
tanLine = derivAt0*theta +1;
plot(theta(L), tanLine(L),'k')

% The second moment is
dd = diff(y,2) ./ dtheta^2;
deriv2At0 = dd(n-1)         % 2nd derivative evaluated at 0
varEst = deriv2At0 - derivAt0^2

% Compare our estimates of mean and variance
[m1 v1] = wblstat(.01,3)
[m2 v2] = wblstat(.02,3)
meanTrue = m1 - m2
varTrue = v2 + v1


%%%%%%%%%%%%%%%%
% Try the same thing again with difference of normals.
%
x = [0:dx:20];
%
pdf1 = normpdf(x,8,1);
pdf2 = normpdf(x,8.5,1);

figure(1), clf
plot(x,pdf1,x,pdf2)
axprefs(gca)
xlabel('x'); ylabel('Probability density')

theta = [-1:dtheta:1]';
% MGF for the difference of the 2nd RV minus the 1st. Notice the
% negative theta on the 1st mgf.
y = mgf(-theta,x,pdf1) .* mgf(theta,x,pdf2) ;

figure(2), clf, hold on
plot(theta,y,'k-','LineWidth',2)
axprefs(gca); xlabel('\theta'); ylabel('M(\theta)')
set(gca,'YLim',[min([.5 .5*min(y)]) 2])
hl = line(get(gca,'XLim'),[1 1])
set(hl,'LineStyle','--','Color','k')


% Use the MGF to get the mean difference. 
d = diff(y) ./ diff(theta);
n = find(theta==0);
% interpolate to approximate the value of the derivative at 0.
derivAt0 = mean([d(n-1) d(n)]);
meanEst = derivAt0
% show the deriv on the fig
L = abs(theta) < .2;
tanLine = derivAt0*theta +1;
plot(theta(L), tanLine(L),'k')

% The second moment is
dd = diff(y,2) ./ dtheta^2;
deriv2At0 = dd(n-1)         % 2nd derivative evaluated at 0
varEst = deriv2At0 - derivAt0^2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3. Relationship between theta_1 and intensity for several "difference"
% distributions. The basic plan here is to consider a pair of random values
% drawn from two distributions. The expectation (and variance) of the RVs
% changes as a function of intensity. The canonical example is that a change in
% intensity leads to a proporational increase and decrease in the two RVs,
% respectively. From these two RVs, we create a difference variable. This
% difference variable has an expetations (trivial) and it has a moment
% generating function. We're interested in the theta_1 term and how it varies as
% a function of the expected difference. 
%
% The reason we care about this is as follows. If the difference RV furnishes
% the increments in a random walk or diffusion, such that the cumulative sum
% begins at 0 and accumulates  until it 
% hits an upper  barrier at +A or a lower barrier at -A, then the probability of
% hitting the A barrier is a logistic function of theta_1 * A.  
%
% If we know how intensity translates to a difference signal, and if we can say
% something about the probability distribution, then we can relate intensity to
% the psychometric function. 
%
% This section gives us a feel for the kind of conditions in which the
% probabilty of stopping at +A as a function of intensity is logistic.

% The following blocks of code illustrate a variety of difference variables that
% might be related to intensity.

% Difference Poissons
theta = [-5:dtheta:5]';     % the axis has to be symmetric about 0.
u1list = [10:.5:14]';
u2list = [10:-.5:6]';
meandiffList = u1list - u2list
figure(1), clf, hold on
for i = 1:length(u1list)
    pdf1 = poisspdf(x,u1list(i));
    pdf1 = pdf1 / trapz(x,pdf1);
    pdf2 = poisspdf(x,u2list(i));
    pdf2 = pdf2 / trapz(x,pdf2);
    
    MGF1 = mgf(theta,x,pdf1);
    MGF2 = mgf(theta,x,pdf2);
    MGFd = MGF1 .* flipud(MGF2);
    plot(theta,MGFd);
    set(gca, 'YLim', [.5 2])
    % Get the root. We're dealing with positive means. So we know where to look
    if meandiffList(i) == 0
        theta1(i) = 0;
    else
        L = theta < -.1;
        k = find(abs((MGFd.*L)-1)  == min(abs((MGFd.*L)-1)))
        theta1(i) = theta(k);
        line([theta(k) theta(k)], [.5 2])
        line(get(gca,'XLim'), [1 1])
    end
    drawnow
end
figure(2),clf, hold on
plot(u1list-u2list,-theta1,'o')
lsline
axprefs(gca)
xlabel('Mean difference')
ylabel('\theta_1')
title('\Delta Poisson: means scale with intensity')

% Difference of gaussians. Assume variance scales w/ mean
theta = [-5:dtheta:5]';     % the axis has to be symmetric about 0.
u1list = [10:.5:14]';
u2list = [10:-.5:6]';
fano = .3
s1list = sqrt(fano * u1list)
s2list = sqrt(fano * u2list)

meandiffList = u1list - u2list
figure(1), clf, hold on
for i = 1:length(u1list)
    pdf1 = normpdf(x,u1list(i),s1list(i));
    pdf1 = pdf1 / trapz(x,pdf1);
    pdf2 = normpdf(x,u2list(i),s2list(i));
    pdf2 = pdf2 / trapz(x,pdf2);
    
    MGF1 = mgf(theta,x,pdf1);
    MGF2 = mgf(theta,x,pdf2);
    MGFd = MGF1 .* flipud(MGF2);
    plot(theta,MGFd);
    set(gca, 'YLim', [0 2])
    % Get the root. We're dealing with positive means. So we know where to look
    if meandiffList(i) == 0
        theta1(i) = 0;
    else
        L = theta < -.1;
        k = find(abs((MGFd.*L)-1)  == min(abs((MGFd.*L)-1)))
        theta1(i) = theta(k);
        line([theta(k) theta(k)], get(gca,'YLim'))
        line(get(gca,'XLim'), [1 1])
    end
    drawnow
end
figure(2),clf, hold on
plot(u1list-u2list,-theta1,'o')
lsline
axprefs(gca)
xlabel('Mean difference')
ylabel('\theta_1')
title('\Delta Gaussian: means scale symmetrically with intensity, var pop to mean')


% Difference of gaussians, under the assumption of assymetric changes to the positively and
% negatively changing means (Pref and Null directions). Assume variance scales
% w/ mean. This is an important case, because it represents a departure from
% conditions that make theta_1 scale with intensity. Nevertheless, it's awfully
% close.
theta = [-5:dtheta:5]';     % the axis has to be symmetric about 0.
u1list = [10:.5:14]'
nfact = .3   % Intensity affects the declining mean by nfact relative the the increasing mean. Use nfact=1 for symmetric changes. 
u2list = u1list(1) + nfact * (u1list(1) - u1list)
fano = .3
s1list = sqrt(fano * u1list)
s2list = sqrt(fano * u2list)

meandiffList = u1list - u2list
figure(1), clf, hold on
for i = 1:length(u1list)
    pdf1 = normpdf(x,u1list(i),s1list(i));
    pdf1 = pdf1 / trapz(x,pdf1);
    pdf2 = normpdf(x,u2list(i),s2list(i));
    pdf2 = pdf2 / trapz(x,pdf2);
    
    MGF1 = mgf(theta,x,pdf1);
    MGF2 = mgf(theta,x,pdf2);
    MGFd = MGF1 .* flipud(MGF2);
    plot(theta,MGFd);
    set(gca, 'YLim', [0 2])
    % Get the root. We're dealing with positive means. So we know where to look
    if meandiffList(i) == 0
        theta1(i) = 0;
    else
        L = theta < -.1;
        k = find(abs((MGFd.*L)-1)  == min(abs((MGFd.*L)-1)))
        theta1(i) = theta(k);
        line([theta(k) theta(k)], get(gca,'YLim'))
        line(get(gca,'XLim'), [1 1])
    end
    drawnow
end
figure(2),clf, hold on
plot(u1list-u2list,-theta1,'o')
lsline
axprefs(gca)
xlabel('Mean difference')
ylabel('\theta_1')
title('\Delta Gaussian: means scale asymmetrically with intensity, var pop to mean')



% Difference of lognormal distributions. I chose lognormals to illustrate an
% extreme case.  Look at the pdfs. I have yet to do the algebra to see why
% theta_1 is linear in the difference variable, but I suspect it's trivial.
% Sorry to be so irresponsible. 

theta = .1 * [-5:dtheta:5]';     % the axis has to be symmetric about 0.
u1list = [2: .05: 2.4]'
nfact = 1   % Intensity affects the declining mean by nfact relative the the increasing mean. Use nfact=1 for symmetric changes. 
u2list = u1list(1) + nfact * (u1list(1) - u1list)
meandiffList = u1list - u2list
theta1 = repmat(nan, size(meandiffList));

figure(1), clf, hold on
for i = 1:length(u1list)

    pdf1 = lognpdf(x,u1list(i),1);
    pdf1 = pdf1 / trapz(x,pdf1);
    pdf2 = lognpdf(x,u2list(i),1);
    pdf2 = pdf2 / trapz(x,pdf2);
    figure(3), clf, hold on
    plot(x,pdf1, x,pdf2)
    
    MGF1 = mgf(theta,x,pdf1);
    MGF2 = mgf(theta,x,pdf2);
    MGFd = MGF1 .* flipud(MGF2);
    figure(1)
    plot(theta,MGFd);
    set(gca, 'YLim', [0 2])
    % Get the root. We're dealing with positive means. So we know where to look
    if meandiffList(i) == 0
        theta1(i) = 0;
    else
        L = theta < -.01;
        k = find(abs((MGFd.*L)-1)  == min(abs((MGFd.*L)-1)))
        theta1(i) = theta(k);
        line([theta(k) theta(k)], get(gca,'YLim'))
        line(get(gca,'XLim'), [1 1])
    end
    drawnow
end
figure(2),clf, hold on
plot(meandiffList,-theta1,'o')
lsline
axprefs(gca)
xlabel('Mean difference')
ylabel('\theta_1')
title('\Delta log: means scale asymmetrically with intensity, var pop to mean')


% Difference between uniform distributions. We'll leave the subtracted
% distribution alone.
theta =  10*[-5:dtheta:5]';     % the axis has to be symmetric about 0.
b1list  = [2:6]';
a1list = b1list - 5;
a2list = repmat(a1list(1),size(a1list));
b2list = repmat(b1list(1),size(b1list));

meandiffList = (b1list+a1list)/2 - (b2list+a2list)/2
theta1 = repmat(nan, size(meandiffList));

figure(1), clf, hold on
for i = 1:length(a1list)

    pdf1 = unifpdf(x,a1list(i),b1list(i));
    pdf1 = pdf1 / trapz(x,pdf1);
    pdf2 = unifpdf(x,a2list(i),b2list(i));
    pdf2 = pdf2 / trapz(x,pdf2);
    figure(3), clf, hold on
    plot(x,pdf1, x,pdf2)
    
    MGF1 = mgf(theta,x,pdf1);
    MGF2 = mgf(theta,x,pdf2);
    MGFd = MGF1 .* flipud(MGF2);
    figure(1)
    plot(theta,MGFd);
    set(gca, 'YLim', [0 2])
    % Get the root. We're dealing with positive means. So we know where to look
    if meandiffList(i) == 0
        theta1(i) = 0;
    else
        L = theta < -.01;
        k = find(abs((MGFd.*L)-1)  == min(abs((MGFd.*L)-1)))
        theta1(i) = theta(k);
        line([theta(k) theta(k)], get(gca,'YLim'))
        line(get(gca,'XLim'), [1 1])
    end
    drawnow
end
figure(2),clf, hold on
plot(meandiffList,-theta1,'o')
lsline
axprefs(gca)
xlabel('Mean difference')
ylabel('\theta_1')


% Here's a bizarre case. Another difference in two uniform distributoins. But all
% we'll do is change the upper bound of one of the distributinos used to form
% the difference. This breaks the theta_1 vs. intensity relationship.
theta =  [-5:dtheta:5]';     % the axis has to be symmetric about 0.
b1list  = [2:6]';
a1list = repmat(0,size(b1list));
a2list = a1list;
b2list = repmat(b1list(1),size(b1list));

meandiffList = (b1list+a1list)/2 - (b2list+a2list)/2
theta1 = repmat(nan, size(meandiffList));

figure(1), clf, hold on
for i = 1:length(a1list)

    pdf1 = unifpdf(x,a1list(i),b1list(i));
    pdf1 = pdf1 / trapz(x,pdf1);
    pdf2 = unifpdf(x,a2list(i),b2list(i));
    pdf2 = pdf2 / trapz(x,pdf2);
    figure(3), clf, hold on
    plot(x,pdf1, x,pdf2)
    
    MGF1 = mgf(theta,x,pdf1);
    MGF2 = mgf(theta,x,pdf2);
    MGFd = MGF1 .* flipud(MGF2);
    figure(1)
    plot(theta,MGFd);
    set(gca, 'YLim', [0 2])
    % Get the root. We're dealing with positive means. So we know where to look
    if meandiffList(i) == 0
        theta1(i) = 0;
    else
        L = theta < -.01;
        k = find(abs((MGFd.*L)-1)  == min(abs((MGFd.*L)-1)))
        theta1(i) = theta(k);
        line([theta(k) theta(k)], get(gca,'YLim'))
        line(get(gca,'XLim'), [1 1])
    end
    drawnow
end
figure(2),clf, hold on
plot(meandiffList,-theta1,'o')
lsline
axprefs(gca)
xlabel('Mean difference')
ylabel('\theta_1')

