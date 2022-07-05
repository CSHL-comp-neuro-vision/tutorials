%% NormalizationAttention_Tutorial.m

addpath('dependencies_NormalizationAttention');
clear all; close all;

%%
% This tutorial works through the 'normalization model of attention' based
% on the 2009 paper published in Neuron by John Reynolds and David Heeger.
% The model is worth going through because it hits on a range of topics
% associated with this course - linear filtering, normalization, population
% responses, spatial attention and feature-based attention.  It's powerful
% because it can predict a range of physiological, behavioral and fMRI
% results from attention studies.
%
% It requires adding the folder 'NormalizationAttention' to the path
%
% Written by G.M. Boynton in June 2010 based heavily on code provided by the
% authors.

%% The 'neural image'

% The backbone of the model is the way the population of neuronal responses
% are represented.  The 'neural image' is a 2-D matrix (or image) in which
% the two matrix dimensions represent two physical dimensions - space and
% orientation in this model - with the intensity of each pixel representing
% the size of the response to a neuron with the corresponding space and
% feature selectivity.   
%
% Showing this matrix as a gray-scale image is a convenient way of
% representing the population response for neurons tuned to specific
% locations (receptive fields) and orientations.

%% Model/neural parameters

% First we'll define our 2-d stimulus space for the neural images.  Space
% (matrix rows) occupies a single dimension which may seem weird, but
% without loss of any generality we can think of the space dimesion as 2-d
% space unwrapped into a single vector.  Or, we can just think of it as 1-d
% space, since all computations will generalize to 2-d space.  If you
% really wanted, you could define a cube, or a 'neural volume' instead of a
% neural image with 2 dimension for space and the third dimension as the
% feature. (and so on).
%
% The structure 'p' will eventually contain all parameters of the model.

clear p
p.x = -200:4:200;              %list of stimulus positions (arbitrary coordinates)
p.theta = -180:1:180;          %list of orientations (degrees)

% The values '4' and '1' determine the step size of the sampling across space and
% orientation.  '4' is a bit coarse but it allows for much faster computing
% time than a more accurate step size of 1.

%% Stimulus parameters

% Stimuli will be defined as having centers and Gassian widths in both
% space and orientation.  These are basically like Gabors if the
% orientation width is narrow and if we assume that the stimulus is
% narrowly defined in spatial frequency.

% We'll define the stimulus parameters inside the structure 'stim'.  Each
% of the fields contain vectors of equal length, with each component
% corresponding to a different stimulus.  We'll start with a single
% stimulus - a high-contrast Gabor at position x = -100 with an orientation
% of theta = 0.

stim.x.center = -100;     %stimulus center positions
stim.x.width = 3;           %stimulus widths (deg)

stim.theta.center = 0;      %stimulus orientation centers
stim.theta.width = 1;       %stimulus orientation widths

stim.contrast = 1;

% A stimulus can be represented as a 'stimulus image' in the same
% coordinates as the neural image.  Each pixel in the stimulus image
% represents the amount of contrast at that pixel's corresponding position
% and orientation.

% We can translate these stimulus parameters into a stimulus image with the
% following function:

img = makeNeuralImage(p,stim);

% You can open the function to see how it works.  It generates Gaussians in
% space and orientation for each stimulus component and computes the outer
% product of the two to generate the stimulus image matrix.  So the
% contrast energy of the stimulus is defined to be a separable matrix with
% Gaussian orientation and spatial tuning on the marginals.

% The stimlulus image (and subsequent neural images) can be visualized with
% the function 'showIt':

figure(1)
clf
showIt(p,img);

% Take some time to think about what different stimuli look like in this
% stimulus space. What does a full-field grating look like?  It's a
% horizontal stripe.  A spatially localized patch of noise?  A vertical
% stripe.  (Note that a 'stimulus image' here is different from an actual
% image of the stimulus)

%% Excitatory response parameters (neuronal tuning widths)
%
% The first stage of the model is the linear excitatory response to the
% stimulus.  This is simply the response of a linear receptive field tuned
% to space and feature.  The center of the tuning for each neuron is
% determined by where it lives in the neural image.  The width of tuning
% for space (receptive field size) and feature (orientation tuning) are
% defined here.  These are in terms of standard deviations of Gaussians.

p.e.x.width = 5;                %RF width
p.e.theta.width = 60;           %orientation tuning width (60)

%% The Excitatory response
%
% The neural response for the excitatory component is a linear filtering of
% the stimulus with the excitatory receptive fields. For any given pixel in
% the neural image, the response is a Gaussian centered at that position
% with widths determined by p.e above, multiplied by the stimulus image,
% and added up.  The entire neural image is therefore calculated with a
% convolution of the stimulus image by the Gaussians determined by the
% excitatory parameters.  

% I've provided a function 'convolveImage' that produces this neural image
% through convolution:

E = convolveImage(p,img,p.e.x.width,p.e.theta.width);

% Here's what the excitatory neural image looks like (the authors call it
% 'excitatory drive'

showIt(p,E);
title('Excitatory drive');

% Does this make sense?  Remember, the intensity of the image at each point
% is the response of a neuron with tuning centered at the corresponding
% location and orientation.  The orientation tuning with is pretty broad
% (60 deg), so neurons tuned fairly far away from the stimulus still
% respond somewhat to the stimulus.  But the spatial receptive field is
% narrow (5 deg), so the neural image drops off rapidly in the spatial (x)
% dimension.

%% Inhibitory response parameters
% 
% In the standard normalization model, the excitatory signal is divided by
% the pooled response across a range of neurons tuned across space and
% features.  We can define the range of spatial pooling the same way we
% defined the stimulus image and the excitatory response parameters:

% Inhibitory (pooling) parameters
p.i.x.width = 20;               %spatial pooling width
p.i.theta.width =360;           %orientation pooling width

% We also need a constant in the denominator to keep the ratio from blowing
% up:
p.sigma = 1e-6;                 %sigma in denominator for normalization

% These specific parameters mean that each neuron is suppressed by neurons
% tuned away across a fairly narrow range in space, but across neurons
% tuned across all orientations.  
%
% The neural image describing the summed response in the pool of neurons is
% therefore simply the convolution of the excitatory neural image (E) with
% Gaussians determined by the inhibitory pooling parameters.

I = convolveImage(p,E,p.i.x.width,p.i.theta.width);

% Here's a picture of the inhibitory pooling neural image which the authors
% call 'inhibitory drive'

showIt(p,I)
title('Inhibitory drive');

% The brightness of each pixel in this image represents the strength of
% divisive inhibition for the corresponding neuron.
%
% Although these inhibitory parameters are defined the same way as the
% excitatory parameters, they have a very different meaning.  The
% excitatory parameters determine the receptive field properties of the
% neurons.  These inhibitory parameters determine the range of neurons that
% we're polling across for normalization. It's similiar in that the
% excitatory parameters determine the range of stimuli excite a given
% neuron, the inhibitory parameters determine the range of neurons that
% inhibit a given neuron.

%%

% The neural image for the normalized response is calculated as:
R = E./(I+p.sigma);

figure(4)
showIt(p,R);
title('Population Response');

% This is the population response for the normalization model.  As
% expected, there are the largest responses for neurons tuned for the
% stimulus properties, and the neural responses fall of for neurons tuned
% away for both space and orientation.  
%
% You should play with the model parameters to see how they affect is
% neural image.
%
% Before we get in to attention, we can play around a bit with this basic
% normalization model to see how it predicts some basic response properties
% of V1 neurons.

%% Contrast response

%Let's measure the response of the normalization model to a range of
%contrasts.  This means generating a neural image for each contrast. I've
%made a function 'normalizationModel' that calculates the neural image, R,
%shown above (and returns the other images for free)

contrastList = exp(linspace(log(1e-5),log(1),21));

R = zeros(length(p.theta),length(p.x),length(contrastList));

for i=1:length(contrastList)
    stim.contrast = contrastList(i);
    R(:,:,i) = normalizationModel(p,stim);
end

%%
% Most of the neurons in these neural images are not responding because
% they are tuned away from the stimulus in space or orientation. Let's find
% the neuron that is most closely tuned to the center of the stimulus and
% plot its contrast response function.

[tmp,xid] = min( (p.x-stim.x.center).^2);
[tmp,thetaid] = min( (p.theta-stim.theta.center).^2);

% and plot the contrast response function:
y1 = squeeze(R(thetaid,xid,:));

figure(2) 
clf
plot(log(100*contrastList),y1,'b-','LineWidth',2);
logx2raw
xlabel('Contrast (%)')
ylabel('Response')
set(gca,'YLim',[0,max(y1)*1.1])

%% Cross-orientation inhibition

% Consider the contrast response for the same neuron, but in the presence
% of a high contrast grating tuned to the orthognal dimension.  This is a
% sequence of plaids where the off-orientation is high contrast and the
% preferred orientation varies:

stim.x.center = [-100,-100];
stim.x.width = [3,3];

stim.theta.center = [0,90];
stim.theta.width = [1,1];

for i=1:length(contrastList)
    stim.contrast = [contrastList(i),1];
    R(:,:,i) = normalizationModel(p,stim);
end

% and plot the contrast response function:

y2 = squeeze(R(thetaid,xid,:));

figure(2)
clf

plot(log(100*contrastList),[y1,y2],'LineWidth',2);
logx2raw
xlabel('Contrast (%)')
ylabel('Response')
set(gca,'YLim',[0,max(y1)*1.1])
legend({'grating','plaid'},'Location','NorthWest');

% Notice how the presence of a high contrast orthogonal stimulus (green)
% suppresses the responses to the stimulus compared to the grating alone
% (blue).  This is exactly what's found in the physiology literature and
% has been modeled by normalization models such as Heeger's work in the
% 90's. Technically, Heeger's original model squares the excitatory input
% before feeding into the normalization process, but the general idea is
% the same.
%
% The normalization model predicts these curves because the normalization
% pool is orientation-independent, so the orthogonal grating adds a strong
% divisive signal.
%
% You probably noticed that the response to the plaid is greater than the
% response to the grating at low grating contrasts.  This is because the
% excitatory orientation tuning is broad (60 deg), so the orthogonal
% stimulus feeds some signal into the excitatory input.

%% Stimulus parameters for an attention experiment

% Finally we're ready to talk about attention.  First we'll set up a new
% stimulus condition.  This will be a classic spatial attention condition
% where two intermediate contrast Gabors are presented, one on the left and
% one on the right side of the visual field.

stim.x.center = [-100,100];
stim.x.width = [3,3];

stim.theta.center = [0,0];
stim.theta.width = [1,1];

stim.contrast = [.25,.25];

%Here's the stimulus image

S = makeNeuralImage(p,stim);
figure(1)
showIt(p,S);
title('Stimulus Image');

%% Attention parameters

% The normalization model of attention defines the spatial and featural
% spread of attention using the same sort of structure as the stimulus,
% excitatory and inhibitory parameters.  

% Spatial attention is a Gaussian 'spotlight' centered at some location
% with some standard deviation that may vary with the task. We'll have 
% attention directed to the left stimulus (x = 100) and focused down to a
% size that matches the size of the stimulus (width = 3) 

attend.x.center = -100;        %center of spatial focus of attention
attend.x.width = 3;            %width of spatial focus of attention 

% Feature-based attention is defined similarly with a center and a width.
% Here's an example of feature-based attention being spread across all
% orientations (width = inf).

attend.theta.center = 0;        %center of attention to orientation
attend.theta.width = inf;       %width of attention to orientation

% The final attention parameters determine the minimum and maximum gain
% changes.  
attend.range = [1,2];

% The 'spotlight' of attention can be represented as an 'Attention Field'
% in the neural image space and can be generated with 'makeNeuralImage'

A = makeNeuralImage(p,attend);
figure(1)
clf
showIt(p,A);
title('Attention Field');

%% Modelling attention by modulating the excitatory input

% The effects of attention are implemented by simply multiplying the
% excitatory drive by the attentional spotlight. 

% Here's the excitatory input, as before:
E = convolveImage(p,S,p.e.x.width,p.e.theta.width);

% And here's the gain change due to attention:
G = A.*E;

figure(1)
showIt(p,G);
title('Excitatory drive with attention');

% See how the excitatory response to the left stimulus is greater
% (brighter) than the right, because attention was directed there.

%% Inhibitory image with attention

% Remakably, the rest of the model is exactly as before.  We calculate the
% inhibitory neural image by convolving the (now modulated) excitatory
% input:

I = convolveImage(p,G,p.i.x.width,p.i.theta.width);

showIt(p,I)
title('Inhibitory drive with attention');

% Notice that for this example, the divisive inhibitory input is also
% greater for the attended stimulus. We'll see soon how the relative
% contributions of attention to the excitatory and inhibitory inputs allows
% for a range of ways that attention can influence the neuronal response.

%% Normalization model with attention

%The last step is also like before, we divide the response by the
%inhibitory input (plus a small constant).
R = G./(I+p.sigma);

showIt(p,R);
title('Population response');

%% Contrast response for attended and unattended stimuli

% The way attention influences neuronal responses as a function of stimulus
% contrast is a useful way to characterize the overall effects of
% attention.

% Next we'll plot contrast response functions for neurons that were
% presented identical physical stimuli, but with attention directed within
% only one of the two receptive fields.  The responses of these two neurons
% is equivalent to the response of a single neuron with attention shifted
% within and away from its receptive field.  

% To plot the contrast response functions, we need to find the indices for
% the two neurons that are most closely tuned to the two stimuli:

for i=1:length(stim.x.center)
    [foo,id] = min( (p.x-stim.x.center(i)).^2);
    xid(i) = id;
    [foo,id] = min( (p.theta-stim.theta.center(i)).^2);
    thetaid(i) = id;
end

% Now we can calculate the neural images and plot the two responses just
% like before.  The function 'normalizationModel' can take in a third
% argument 'attend' that contains the attention parameters.  The
% calculations within 'normalizationModel' are identical to those in this
% script.

R = zeros(length(p.theta),length(p.x),length(contrastList));

for i=1:length(contrastList)
    stim.contrast = contrastList(i)*[1,1];
    R(:,:,i) = normalizationModel(p,stim,attend);
end

% We can pull out the two contrast response functions to each stimulus from
% the neural images:

y = zeros(length(contrastList),length(xid));
for i=1:length(xid)
    y(:,i) = R(thetaid(i),xid(i),:);
end

figure(2)
clf
plot(log(100*contrastList),y,'LineWidth',2);
logx2raw
xlabel('Contrast (%)')
ylabel('Response');

legend({'Attend in','Attend out'},'Location','NorthWest')

% There it is, a higher response to an attended stimulus than an unattended
% stimulus.  Note that for these parameters, it looks like attention is
% acting as a 'response gain', which is a vertical scaling between the
% attended and unattended responses across contrast.  

%% Response gain vs. Contrast gain.

% Why does a narrow focus of attention lead to response gain?  It helps to
% look at the cross-section of the neural images to get some insight.
% Consider a high-contrast stimulus with a narrow focus of attention. We'll
% make the conditions more extreme (larger stimulus, smaller focus of
% attention) for better illustration.  Here is a complete set of stimulus
% and attentional parameters:

stim.x.center = [-100,100];
stim.x.width = [20,20];           %stimulus widths (deg)

stim.theta.center = [0 0];
stim.theta.width = [1,1];

stim.contrast= [1,1];

attend.x.center = -100;
attend.x.width = 1;         %narrow!

attend.theta.center = 0;
attend.theta.width = inf;

attend.range = [1,2];

[R,S,E,A,G,I] = normalizationModel(p,stim,attend);
%(The function can return all of the neural images involved in the
%computation)

%Here's a plot of the spatial profile of the exciatory drive and inhibitory
%drives after it is modulated by attention (G).

id = round(length(p.theta)/2);  %find the middle row of neural image for plotting
figure(2)
clf
subplot(2,1,1)
plot(p.x,G(id,:),'g-','LineWidth',2);

xlabel('Space');
title('Excitatory drive');

subplot(2,1,2)
plot(p.x,I(id,:),'g-','LineWidth',2);
title('Inhibitory drive');
xlabel('Space');

% You can see how the effect of attention on the excitatory drive is to
% multiply the original excitatory input by a narrow Gaussian (ranging
% between 1 and 2), which increases the excitatory drive only within cells
% with RF's near the attended location (like the one we're plotted from
% above).

% The inhibitory drive is the convolution of the excitatory drive by a
% pooling filter, which simply spatially blurrs the excitatory drive. Since
% the spike is so narrow, the spatial blurring by the pooling process
% leaves the inhibitory drive very similar for the attended and unattened
% regions of space (left vs. right).

% The output of the model is basically the ratio of the excitatory and the
% inhibitory drives.  Since convolution is linear, changing the contrast
% simply scales these curves up and down (try running this section again
% with a different contrast.  Only the y-axis scales). So you can see how
% for the neurons with RF's centered at the focus of attention, the effect
% of attention is all in the numerator.  So scaling by contrast simply
% scales the response.  Contrast gain!

%%

% Now, consider the same stimulus but with a broad focus of spatial attention:
attend.x.width = 30;

% and overlay the spatial profile of excitatory and inhibitory drives:
[R,S,E,A,G,I] = normalizationModel(p,stim,attend);

figure(2)
subplot(2,1,1)
hold on
plot(p.x,G(id,:),'b:','LineWidth',2);

subplot(2,1,2)
hold on
plot(p.x,I(id,:),'b:','LineWidth',2);
legend({'Narrow focus','Broad focus'});

% The effect of a broad focus of spatial attention is to boost a broader 
% neuronal population since we multiplied by a broader Gaussian.  Note,
% however, that the peak excitatory drive is the same as for the narrow
% focus. 

% This time, the inhibitory drive is strongly affected by attention.  This
% is because the blurring of the new excitatory drive by the attention
% filter is summing the response over a lot of active neurons. 

% As before, contrast simply scales both the excitatory and inhibitory
% drives up and down.  But now both the numerator AND denominator for the
% model are growing faster with contrast for the 'attended' neuron.  This
% is just like changing the contrast of the stimulus with attention.
% Contrast gain!  

% Reynolds and Heeger go and predict a range of published papers showing
% both contrast gain and response gain effects of attention just by
% changing the spatial focus of attention relative to the size of the
% stimulus. This helps to settle a lot of discrepancies (and debates) in
% the literature.  This simple explanation also leads to some testable
% hypotheses about how the spatial focus of attention should affect
% neronal, behavioral and fMRI measurements.

%%  Attention and orientation tuning.

% McAdams and Maunsell (1999) measured how the orientation tuning of V4
% neurons are affected by spatial attention by measuring orientation tuning
% while monkeys attended either inside or outside the receptive field of
% each neuron. 

% The model naturally predicts how spatial and feature-based attention
% influences orientation tuning curves.  We can see this by simply slicing
% through the neural images along the orientation dimension.  For
% simplicity, we'll assume that as before, while feature-based attention is
% localized, attention is directed to all orientations.

% Let's set up some reasonable stimulus and attention conditions:

stim.x.center = [-100,100];
stim.x.width = [10,10];           %stimulus widths (deg)

stim.theta.center = [0 0];
stim.theta.width = [1,1];

stim.contrast= [1,1];

attend.x.center = -100;
attend.x.width = 10;

attend.theta.center = 0;
attend.theta.width = inf;

attend.range = [1,4];

R = normalizationModel(p,stim,attend);

figure(2)
clf
plot(p.theta,R(:,xid(1)),'r-','LineWidth',2)
hold on
plot(p.theta,R(:,xid(2)),'b-','LineWidth',2);
xlabel('Orientation (deg)')
legend({'attended','unattended'});

% This looks just like the scaling of the tuning functions seen by McAdams and Maunsell
% (1999).  

% You might think it's weird that we're plotting slices of the neural image
% which is really the model's prediction of the response of a bunch of
% different neurons to the same stimulus.  But if you think about it, this
% is just the same as plotting the response of the same neuron to a range
% of stimuli.  



% The model (as I've implemented it) used this 'or' rule all along but it
% wasn't noticeable because attend.theta.width = inf, so 'and' and 'or'
% predict the same thing.

% But with attend.theta.width = 30, we can now look at the effects of
% feature based attention.

%% Moran and Desimone

% The study by Moran and Desimone (1985) was the first to show attentional
% effects in monkey cortex.  V4 neurons were measured with two oriented bars 
% in the receptive field.

stim.x.center = [-100,-100];
stim.x.width = [10,10];           %stimulus widths (deg)

stim.theta.center = [0 90];
stim.theta.width = [1,1];

stim.contrast= [1,1];

attend.x.center = -100;  %attend left
attend.x.width = 2;

attend.theta.center = 0; %attend horizontal
attend.theta.width = 30;

attend.range = [1,2];
attend.method = 'add';

Rpref  = normalizationModel(p,stim,attend);

attend.theta.center = 90;  %attend vertical
Rnonpref = normalizationModel(p,stim,attend);

y1= Rpref(thetaid(1),xid(1));
y2 =Rnonpref(thetaid(1),xid(1));

figure(2)
clf
bar([y1,y2]);
set(gca,'XTickLabel',{'attend pref','attend non-pref'});
set(gca,'XLim',[.5,2.5])
ylabel('Model Response');

%% Feature-similarity gain

% This 'cross' profile of attentional gain means that feature-based
% attention should influence responses for neurons with receptive fields
% outside the focus of spatial attention.  This turns out to be true.
% Treue and Martinez-Trujillo (1999) originally discovered this by
% measuring the response to an unattended moving stimulus in MT when
% feature-based attention was altered by having the monkey perform a
% task on a stimulus in the opposite hemfield.  

% We can demonstrate this with the model by presenting stimuli that simulate
% the conditions for an fMRI experiment by Saenz et al. (2002).

stim.x.center = [-100,-100,100];
stim.x.width = [10,10,10];

stim.theta.center = [0,90,0];
stim.theta.width = [2,2,2];
stim.contrast = [1,1,1];

attend.x.center = -100;
attend.x.width = 2;
attend.theta.center = 0;
attend.theta.width = 30;

[Rpref,S,E,Apref,G,I] = normalizationModel(p,stim,attend);

% Let's look at the stimulus:
figure(1)
showIt(p,S);
title('Stimulus image');

% In the real experiment, the stimuli were moving dots instead of oriented
% stimuli but it doesn't matter here.  On the attended (left) side, two
% stimuli were presented overlapping in space - one vertical and one
% horizontal.  The unattended side contained a single oriented stimulus.

% Save the response to the unattended stimulus for the neuron that both
% prefers that stimulus and the attended orientation (0).

yP = Rpref(thetaid(2),xid(2));

% Now we'll shift feature-based attention to the non-preferred component of
% the attended stimulus.

attend.theta.center = 90;


Rnonpref = normalizationModel(p,stim,attend);

% Pull out the same neuron's response
yNP = Rnonpref(thetaid(2),xid(2));

% And plot both
figure(2)
clf
bar([yP,yNP])
set(gca,'XLim',[.5,2.5]);
set(gca,'XTickLabel',{'Attend out, preferred','Attend out, non-preferred'});
ylabel('Model response');

