%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shTutorial1.m - Introduction and overview to the MT model
%% Authors: Timothy Saint and Eero Simoncelli
%%
%% This file is part of the MTmodel package, available at
%%      http://www.cns.nyu.edu/~lcv/MT-model.html
%% See the README file for brief description and installation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This tutorial provides a basic introduction to a package of matlab
% code implementing a two-stage functional model for neurons in visual
% areas V1 and MT, as published in "A Model of Neuronal Responses in
% Visual Area MT", EP Simoncelli and DJ Heeger, Vision Research,
% 38(5):743-761, March 1998.  The tutorial describes visual stimulus
% generation, parameters and options, computation of model V1 or MT
% responses, and generation of basic tuning curves.  It concludes with
% a sequence of examples that generate figures from the original
% publication.  The model is further described in a second tutorial,
% shTutorial2.m.  

% We recommend that you use the tutorial by viewing it in a text
% editor (we recommend matlab's editor or gnu emacs), reading through
% the commented paragraphs, and executing each chunk of matlab code
% (using copy-paste) as you go.

% ----------------------------------------------------------------------
% I. SETUP
% 
% Follow the installation instructions in the README file.  Make sure
% that the directory containing the MTmodel code is in your path.  We
% recommend:

% addpath(genpath(MYPATH));

% where MYPATH is the location of the MTmodel folder.  If your path is
% correct, you should be able to execute the following to get a quick
% list of functions in the package, along with tutorial files:

help MTmodel

% ----------------------------------------------------------------------
% II. MODEL BASICS: INPUTS, PARAMETERS, OUTPUTS

% The input for the model is a visual stimulus sequence (i.e., a
% movie), represented in matlab as a 3D array.  We provide a number of
% functions for generating standard stimuli, such as drifting
% sinusoidal gratings or dots.  For example:

help mkDots
dir = 0;  speed = 1;
s = mkDots([39 39 50], dir,  speed,  0.15);

% The first argument specifies the size of the movie as [Y X T].  This
% may seem like an odd ordering for the dimensions, but it is
% consistent with matlab's 2D matrix representation [Y X].

% You can view a matlab movie of the stimulus with the FLIPBOOK
% function. You will probably see some temporal irregularities when
% viewing the movie.  These do not indicate problems with the stimulus
% -- they arise because flipbook loops through the movie displaying
% each frame without attempting to lock to the refresh rate of the
% monitor (or any other timing signal).

flipBook(s);

% The behavior of the model is governed by a set of parameters, which
% are bundled up in a single MATLAB structure.  We provide an
% "original" set of parameters that replicate the model as described
% in the original published article:

pars = shPars

% The various fields of this structure correspond to different parameters,
% each of which is used in one or more stages of the model. Each will be
% explained in detail in later tutorials.

% The model is built in two successive stages, corresponding to
% directionally-selective complex cells in visual area V1, and
% pattern-selective cells in visual area MT.  The responses of neurons
% (over time) at either of these stages may be computed using the
% SHMODEL function:

[pop, ind] = shModel(s, pars, 'mtPattern');

% The first and second arguments are the stimulus and the model
% parameters.  The third argument specifies the stage of the model
% whose output you are requesting.  There are many options that allow
% you to look at outputs of intermediate stages of the computation,
% but the two most important stages are:
%
% 'v1Complex'   - directionally-selective V1 complex cells
% 'mtPattern'   - MT pattern cells 


% ----------------------------------------------------------------------
% IV. EXAMINING RESPONSES OF MODEL MT NEURONS

% The POP output of the call to the SHMODEL command contains the
% responses (over time) of all of the neurons in the population.  This
% includes neurons with receptive fields centered at different spatial
% locations, and with different direction and speed preferences.  The
% direction and speed preferences are specified by a field of the PARS
% structure:

pars.mtPopulationVelocities

% This is a 19x2 matrix, in which each row (a 2-vector) specifies the
% preferred direction and speed, in units of radians and pixels/frame,
% of one of the neurons in the MT population.  The responses of these
% neurons are computed every time you run the model.  Notice that the
% second neuron in the population prefers [0, 1], which is motion to
% the right at the speed of one pixel/frame.  That's identical to the
% speed and direction of the motion of the dots stimulus you created
% above.

% The neural responses are embedded in the output matrix POP.
% Specifically, the j-th row of POP contains the temporal responses of
% all neurons that have the velocity tuning indicated by row j of
% pars.mtPopulationVelocities.  This arrangement is useful for
% computations but not easy to examine.  Therefore we provide the
% function SHGETNEURON, which returns the responses of all of the
% neurons with receptive fields centered at a particular spatial
% position:

p = shGetNeuron(pop, ind);

% By default, SHGETNEURON returns the responses of the neurons in the
% population whose receptive fields are closest to the center of the
% stimulus.

% Each row of the output matrix P contains the response (over time) of
% the neuron whose preferred velocity is specified by the
% corresponding row of pars.mtPopulationVelocities.  For example, the
% second row of P contains the responses of the neuron whose preferred
% velocity is [0 1].  Plotting all of the responses in the same
% figure:

figure(1);
[sortedResp, sortedInd] = sort(mean(p,2),1,'descend'); %sort by mean response
plot([0 0],[1 1],'w');                                % dummy point, for legend title
hold on; plot(p(sortedInd,:)'); hold off;
set(gca,'XLim',[1 size(p,2)]);
xlabel('time (frames)');  ylabel('response (arb. units)');
title(sprintf('MT responses to dots drifting at [%.1f, %.1f]', dir,speed));
vels = pars.mtPopulationVelocities;  nvels = 5;
eval(['legend(gca,''pref. velocity''',sprintf(',''[%.1f, %.1f]''', vels(sortedInd(1:nvels),:)'), ')']);

% The responses are probably a bit irregular, because the random dot
% stimulus is a bit irregular. But you should be able to see that one
% neuron is responding much more than the others.  As the legend
% indicates, this is the neuron tuned for motion in the direction [0,
% 1].  The other neurons respond more or less strongly depending on
% their tuning, but none so much as the one tuned to the speed of this
% stimulus.

% ----------------------------------------------------------------------
% VI. EXAMINING RESPONSES OF MODEL V1 NEURONS

% As mentioned above, it's also possible to look at the output of the
% V1 neurons using the SHMODEL command. And as for MT, the model always
% computes the responses of a population of V1 neurons:

pars.v1PopulationDirections

% This time the matrix is 28x2; each row is again a different neuron
% in the population, but this time the columns have a slightly
% different meaning. The first column is still the direction of motion
% of the neuron's preferred grating (in units of radians with 0
% corresponding to rightward).  The second column is the ratio of the
% temporal frequency to the spatial frequency (in pixels/frame) of the
% neuron's preferred grating.  For example, the first neuron in the
% population prefers motion to the left and slightly up, at a speed
% of 1.5681 cycles per/frame. 

% The version of the model implemented here uses a population of V1
% neurons whose peak spatio-temporal frequency values lie on the
% surface of a sphere in the Fourier domain.  For computational
% efficiency, the standard implementation provided does not include
% neurons at multiple scales.  We can look at the population
% graphically by plotting all their preferred directions, spatial
% frequencies, and temporal frequencies in the 3D Fourier domain:

shShowV1PopulationDirectionsDots(pars.v1PopulationDirections);

% Now we generate a drifting sinusoidal grating stimulus, with random
% direction and speed.  We encourage you to rerun all the code in this
% section a few times to see how the model responds to different
% gratings.

dims = shGetDims(pars, 'v1Complex', [1 1 20]);
direction = 2*pi*rand;  
tfSfRatio = 2*rand; 
gratingParams = v12sin([direction, tfSfRatio]);
spatialFrequency = gratingParams(2);
temporalFrequency = gratingParams(3);
s = mkSin(dims, direction, spatialFrequency, temporalFrequency);
s = (s+1)./2;
flipBook(s);

% A few words about the choice of stimulus dimensions (DIMS).  The
% computational backbone of this model is a series of convolutions
% that are performed in multiple stages.  The stimulus is convolved
% with a set of filters in shModelV1Linear to get the first part of
% the V1 responses.  Then these are squared and blurred (by convolving
% with another filter), etc.  For each of these convolutions, the
% output is smaller (both spatially and temporally) than the input.
% So, working backwards, in order to guarantee that you'll end up with
% at least one output response, you have to start with a stimulus of a
% particular size (or larger).  This size depends on the stage of the
% model that you want to compute, and the function SHGETDIMS is
% provided to calculate this size for you. The third argument
% specifies how large an output you want, in [y x t] coordinates.  In
% the example above, we don't care about having the responses at many
% spatial positions, but we do want to see the response over time,
% which is why [1 1 20] is the argument.

% The function GRATINGPARAMS also requires some explanation. As you
% saw above, all the neurons in the population lie on a sphere in
% fourier space. If use as your stimulus a grating with arbitrary
% parameters, chances are it won't lie near the sphere. The function
% V12SIN takes as arguments a direction and a ratio of temporal to
% spatial frequencies, and returns the parameters of a grating that
% lies on the sphere and also has the right direction and TF/SF
% ratio. The function has that name because you can also use it to
% find the preferred direction, SF, and TF of any neuron in the
% population.

% That said, let's look at the actual responses. Again, you'll have to use
% SHGETNEURON to extract the responses of a single set of neurons
% at the same location:

[pop, ind] = shModel(s, pars, 'v1Complex');
p = shGetNeuron(pop, ind)

figure; 
[sortedResp, sortedInd] = sort(mean(p,2),1,'descend');  % sort by mean response
plot([0 0],[1 1],'w');                                  % dummy point, for legend title
hold on; plot(p(sortedInd,:)'); hold off;
set(gca,'XLim',[1 size(p,2)]);
xlabel('time (frames)');  ylabel('response (arb. units)');
title(sprintf('V1 responses to sinusoidal grating drifting at [ %.1f, %.1f]',dir,speed));
dirs = pars.v1PopulationDirections; ndirs= 5;
eval(['legend(gca,''pref. [dir, speed]''',sprintf(',''[%.1f, %.1f]''', dirs(sortedInd(1:ndirs),:)'), ')']);

% The responses are all essentially flat, but the 28 different neurons
% respond at different levels. You should be able to see that the
% neurons with the largest responses are those that have preferred
% [direction speed] closest to that of the stimulus.

% -------------------------------------------------------------------------
% VII. EXAMINING THE RESPONSE OF NEURONS NOT IN THE POPULATION

% The model uses the minimal number of neurons necessary to smoothly
% and completely represent the responses of a dense population.  For
% any one of a number of reasons you might want to compute the
% responses of neurons that aren't in the base (computed) population.
% This can be done in an exact way by interpolating, and this
% capability is included in the implementation of SHMODEL.  You simply
% add a fourth argument:

dir = 0;  speed = 1;
dims = shGetDims(pars, 'mtPattern');
s = mkDots([95 95 31], dir, speed, .15);
tuningVelocities = [0 1; 0 4]; 
[pop, ind, res] = shModel(s, pars, 'mtPattern', tuningVelocities);

% The fourth argument has the same format as
% pars.mtPopulationVelocities: each row specifies a neuron's preferred
% velocity in [direction, speed] format.  In the example above, we've
% asked for neurons that prefer motion to the right, but at different
% speeds (1 and 4 pixels/frame).

% Note that the responses of the full population of neurons are still
% computed, and returned in POP. But now you also get RES, which
% contains the responses of the two neurons you asked for. You can
% look at these responses using the code you used above but
% substituting RES for POP:

r = shGetNeuron(res, ind);

figure
plot([0 0],[1 1],'w');                                % dummy point, for legend title
hold on; plot(r'); hold off
set(gca,'Xlim', [1 size(r,2)]);
title(sprintf('Responses of MT neurons to dots drifting at [%.1f %.1f]',dir,speed));
eval(['legend(gca,''pref. velocity''',sprintf(',''[%.1f, %.1f]''', tuningVelocities'), ')']);

% The neuron with the larger response is the neuron whose tuning
% velocity matches the velocity of the stimulus.  Incidentally, this
% neuron happens to be in the population, which you can see if you
% plot the full population responses (see examples in preceding
% section).

% ----------------------------------------------------------------------
% VII. GENERATING TUNING CURVES
%
% By looping over stimulus parameters and taking the temporal average
% of one cell's response to each stimulus, it's easy to generate
% simulated tuning curves. For example:

v1params = [pi/3, 1];                   % params for a particular v1 neuron
prefGrating = v12sin(v1params);         % calculate the preferred grating of this neuron
dims = shGetDims(pars, 'v1Complex', [1 1 20]);
dirs = linspace(0, 2*pi, 13);           % create vector of stimulus directions to test
r = zeros(size(dirs));
for n = 1:length(dirs)
  fprintf(1,'.');
  s = mkSin(dims, dirs(n), prefGrating(2), prefGrating(3));
  [pop, ind, res] = shModel(s, pars, 'v1Complex', v1params);
  r(n) = mean(shGetNeuron(res, ind));
end
plot(dirs*180/pi, r, 'r-', dirs*180/pi, r, 'k.');
set(gca,'Xlim',[0 360], 'XTick',90*[0:4]);
xlabel('grating direction (deg)');  ylabel('mean response');
ttl = sprintf('grating response, V1 complex cell tuned for %.0fdeg, %.1fpix/frame', v1params(1)*180/pi, v1params(2));
title(ttl);

% The orientation tuning curve indicates that this complex cell is
% highly selective for stimulus direction.  You might prefer to view
% this as a polar plot:

polar(dirs, r, 'r-');
hold on; polar(dirs, r, 'k.'); hold off 
title(ttl);

% We provide a number of functions that can generate tuning curves for
% you automatically, plotting each point as it is computed.  For
% example:

[x, y] = shTuneGratingDirection(pars, [pi/2, 1], 'v1Complex', 12);

% This can be run on any stage of the model just by changing the 3rd
% argument.  A listing of additional tuning curve functions can be
% viewed by executing "help MTmodel", and the source code for these is
% in the "tune" subdirectory.  Here are some more examples:

[x, y] = shTuneGratingSf(pars, [pi/2 1], 'v1Complex', 12);
[x, y] = shTuneDotDirection(pars, [pi/2 1], 'v1Complex', 12);

% ------------------------------------------------------------------------
% VIII. FIGURES FROM THE ORIGINAL JOURNAL ARTICLE

% This last section of the tutorial provides code to generate some of
% the figures from the original article that describes the model (Vision
% Research, 38(5):743-761, 1998).  This is done using commands that
% have been introduced above, or commands that are similar.

pars = shPars;

%--------
figure(9); clf reset; set(gcf, 'color', 'w');

[x, v1Sin] = shTuneGratingDirection(pars, [0 .35], 'v1Complex', 21);
[x2, v1Plaid] = shTunePlaidDirection(pars, [0 .35], 'v1Complex', 21);

[x3, mtSin] = shTuneGratingDirection(pars, [0 .35], 'mtPattern', 21);
[x4, mtPlaid] = shTunePlaidDirection(pars, [0 .35], 'mtPattern', 21);

axisMax = 1.25*max(mtSin);
subplax([2 2], [1.15 1.15], [.7 .7], 1);
polar(0,axisMax); hold on; polar(x, mtSin, 'k.-'); hold off
title('MT, grating') 
%
axisMax = 1.25*max(mtPlaid);
subplax([2 2], [1.15 2.15], [.7 .7], 1);
polar(0,axisMax); hold on; polar(x, mtPlaid, 'k.-'); hold off
title('MT, plaid') 
%
axisMax = 1.25*max(v1Sin);
subplax([2 2], [2.15 1.15], [.7 .7], 1);
polar(0, axisMax); hold on; polar(x, v1Sin, 'k.-'); hold off
title('V1, grating') 
%
axisMax = 1.15*max(v1Plaid);
subplax([2 2], [2.15 2.15], [.7 .7], 1);
polar(0, axisMax); hold on; polar(x, v1Plaid, 'k.-'); hold off
title('V1, plaid') 

%--------
%% Warning: this one is a bit slow...
figure(10); clf reset; set(gcf, 'color', 'w');
currentPos = get(gcf, 'position');
set(gcf, 'position', [currentPos(1), currentPos(2), 500, 900]);

dims = shGetDims(pars, 'mtPattern');

neurons = [0 1.5; 0 .125; 0 9];
speedMinMax = [.3125 5; .0375 .6; 1 10];
barEdgeWidth = [2 1 11];
nDataPoints = 6;
xSpeed = zeros(3, nDataPoints);
yResponseToPreferred = zeros(3, nDataPoints);
yResponseToAntiPreferred = zeros(3, nDataPoints);
nullReponse = zeros(3, nDataPoints);
% Loop over each neuron: first the bandpass neuron, then the lowpass
% neuron, then the highpass neuron.
for iNeuron = 1:3 %3
    % Compute the response of this neuron vs. speed for a bar moving in the
    % preferred direction.
    [xSpeedTmp, yResponseToPreferredTmp] = shTuneBarSpeed(pars, neurons(iNeuron, :), ...
                                  'mtPattern', nDataPoints, speedMinMax(iNeuron, 1), ...
                                  speedMinMax(iNeuron, 2), neurons(iNeuron, 1), ...
                                  1, barEdgeWidth(iNeuron));
    % Compute the response of this neuron vs. speed for a bar moving in the
    % antipreferred direction.
    [xSpeedTmp, yResponseToAntiPreferredTmp] = shTuneBarSpeed(pars, neurons(iNeuron, :), ...
                                  'mtPattern', nDataPoints, speedMinMax(iNeuron, 1), ...
                                  speedMinMax(iNeuron, 2), neurons(iNeuron, 1)+pi, ...
                                  1, barEdgeWidth(iNeuron));                                 
    % Compute the null response of this neuron
    [pop, ind, nullResponseTmp] = shModel(zeros(dims), pars, 'mtPattern', neurons(iNeuron, :));
    nullResponseTmp = mean2(shGetNeuron(nullResponseTmp, ind));
    nullResponseTmp = nullResponseTmp.*ones(1, nDataPoints);
    % Store all the computed responses in matrices that will contain the
    % responses of all three neurons.
    xSpeed(iNeuron, :) = xSpeedTmp;
    yResponseToPreferred(iNeuron, :) = yResponseToPreferredTmp;
    yResponseToAntiPreferred(iNeuron, :) = yResponseToAntiPreferredTmp;
    nullResponse(iNeuron, :) = nullResponseTmp;
end

% Now plot all the responses
subplax([3 1], [3.1 1.1], [.7 .7], 1);
semilogx(xSpeed(1,:), yResponseToPreferred(1,:), 'b-', xSpeed(1,:), yResponseToPreferred(1,:), 'k.');
hold on
semilogx(xSpeed(1,:), yResponseToAntiPreferred(1,:), 'r-', xSpeed(1,:), yResponseToAntiPreferred(1,:), 'k.');
semilogx(xSpeed(1,:), nullResponse(1,:), 'k--');
hold off
title('"bandpass" speed tuning curve');
xlabel('speed (px/frame)');
ylabel('response');
axis([xSpeed(1,1), xSpeed(1,end), 0, 1.2*max(yResponseToPreferred(1,:))]);
%
subplax([3 1], [2.1 1.1], [.7 .7], 1);
semilogx(xSpeed(2,:), yResponseToPreferred(2,:), 'b-', xSpeed(2,:), yResponseToPreferred(2,:), 'k.');
hold on
semilogx(xSpeed(2,:), yResponseToAntiPreferred(2,:), 'r-', xSpeed(2,:), yResponseToAntiPreferred(2,:), 'k.');
semilogx(xSpeed(2,:), nullResponse(2,:), 'k--');
hold off
title('"lowpass" speed tuning curve');
xlabel('speed (px/frame)');
ylabel('response');
axis([xSpeed(2,1), xSpeed(2,end), 0, 1.2*max(yResponseToPreferred(2,:))]);
%
subplax([3 1], [1.1 1.1], [.7 .7], 1);
semilogx(xSpeed(3,:), yResponseToPreferred(3,:), 'b-', xSpeed(3,:), yResponseToPreferred(3,:), 'k.');
hold on
semilogx(xSpeed(3,:), yResponseToAntiPreferred(3,:), 'r-', xSpeed(3,:), yResponseToAntiPreferred(3,:), 'k.');
semilogx(xSpeed(3,:), nullResponse(3,:), 'k--');
hold off
title('"highpass" speed tuning curve');
xlabel('speed (px/frame)');
ylabel('response');
axis([xSpeed(3,1), xSpeed(3,end), 0, 1.2*max(yResponseToPreferred(3,:))]);

%--------
figure(11); clf

[x, yPref] = shTuneDotCoherence(pars, [0 1], 'mtPattern', 8, 71);
[x, yAntipref] = shTuneDotCoherence(pars, [0 1], 'mtPattern', 8, 71, pi);

plot(x, yPref, 'r-', x, yPref, 'k.');
hold on;  plot(x, yAntipref, 'b-', x, yAntipref, 'k.');  hold off
xlabel('dot stimulus coherence');
ylabel('response');

%--------
figure(12); clf

% No prepackaged function this time.  Calculate how big our stimuli
% should be:
dims = shGetDims(pars, 'mtPattern', [1 1 71]);
neuron = [0 1];
pref = neuron;

% Initialize vectors that will store the neuronal responses.
nPref = [0, 8, 16, 64, 256];
nMask = [0, 16, 64, 256];
cols = 'rgbk';
h=zeros(1+length(nMask),2);
yRes = zeros(size(nPref, 2), size(nMask, 2));
rfRad = 15.5;
rfArea = pi.*rfRad.^2;

% Compute the null response of the model neuron.
[pop, ind, resNull] = shModel(zeros(dims), pars, 'mtPattern', neuron);
yNull = mean2(shGetNeuron(resNull, ind));
yNull = yNull.*ones(1, size(yRes, 1));
w = mkWin(dims, 15, 2);
for j = 1:size(nMask, 2)
  for i = 1:size(nPref, 2)
    fprintf(1,'.');
    % Generate the stimuli
    dPref = nPref(i)./rfArea;
    dMask = nMask(j)./rfArea;

    sDots = mkDots(dims, pref(1), pref(2), dPref);
    sMask = mkDots(dims, pref(1)+pi, pref(2), dMask);
    sDotsWithMask = w .* min(sDots + sMask, 1);

    % Compute the response of the model neuron.
    [pop, ind, res] = shModel(10.*sDotsWithMask, pars, 'mtPattern', neuron);
    yRes(i, j) = mean2(shGetNeuron(res, ind));
  end

  % Now display the results so far
  clf;
  x = nPref(1:i);
  h(1,1) = plot([0 0],[1 1],'w');                                % dummy point, for legend title
  hold on
  for n = 1:length(nMask)
    h(n+1,:) = plot(x, yRes(1:i, n), sprintf('%c-',cols(n)), x, yRes(1:i, n), 'k.');
  end
  plot(x, yNull(1:i), 'k--');
  axis([min(x) max(x) 0 1.2*max2(yRes)]);
  hold off
  drawnow
end
title('MT cell response to mixture of preferred and antipreferred');
eval(['legend(h(:,1),''num antipref dots''',sprintf(',''%d''', nMask'), ')']);
xlabel('Number of dots in preferred direction');
ylabel('response');


%--------
figure(13); clf

nFrames = 101;  % nFrames is the number of frames in each stimulus. The higher 
                % nFrames is, the smoother the resulting tuning curve 
                % will be, but the longer the code will take to run.
[x, y] = shTuneDotMaskDirection(pars, [pi 1], 'mtPattern', 9, nFrames);
x = x.*180./pi; % convert from radians to degrees.
x = [x, 360];
y = [y, y(1)];
h1 = plot(x, y, 'r-', x, y, 'k.');
hold on;  
h2 = plot(x, max(y).*ones(size(y)), 'k--');  
hold off
axis([0 360 0 1.2*max(y)]);
title('MT cell response to preferred dots + mask dots');
xlabel('mask direction (degrees)');
ylabel('response');
legend([h1(1),h2(1)], 'preferred + mask', 'preferred alone');


%--------
figure(14); clf

% No prepackaged function this time.  Calculate how big our stimuli
% should be:
dims = shGetDims(pars, 'mtPattern', [1 1 15]);
neuron = [0 1];
pref = neuron;

% Initialize vectors that will store the neuronal responses.
x = linspace(-pi, pi, 9);
yDots = zeros(size(x));
yDotsWithMask = zeros(size(x));

% Compute the null response of the model neuron.
[pop, ind, resNull] = shModel(zeros(dims), pars, 'mtPattern', neuron);
yNull = mean2(shGetNeuron(resNull, ind));
yNull = yNull.*ones(size(yDots));
h = zeros(3,2);
for i = 1:length(x)
  % Generate the stimuli
  sDots = mkDots(dims, x(i), pref(2), .15);
  sMask = mkDots(dims, pref(1)+pi, pref(2), .15);
  sDotsWithMask = sDots + sMask;
  sDotsWithMask(sDotsWithMask > 1) = 1;

  % Compute the response of the model neuron.
  [pop, ind, resDots] = shModel(sDots, pars, 'mtPattern', neuron);
  [pop, ind, resDotsWithMask] = shModel(sDotsWithMask, pars, 'mtPattern', neuron);
  yDots(i) = mean2(shGetNeuron(resDots, ind));
  yDotsWithMask(i) = mean2(shGetNeuron(resDotsWithMask, ind));

  % Now display the results so far
  clf;
  h(1,:) = plot(180*x(1:i)/pi, yDots(1:i), 'r-', 180*x(1:i)/pi, yDots(1:i), 'k.');
  hold on
  h(2,:) = plot(180*x(1:i)/pi, yDotsWithMask(1:i), 'b-', 180*x(1:i)/pi, yDotsWithMask(1:i), 'k.');
  h(3,:) = plot(180*x/pi, yNull, 'k--');
  axis([-180 180 0 1.2*max([yNull, yDots, yDotsWithMask])]);
  hold off
  drawnow
end
title('MT cell, direction tuning for dots');
legend(h(:,1), 'single dot field', 'w/ antipreferred mask', 'spontaneous');
xlabel('Direction (degrees)');
ylabel('response');
