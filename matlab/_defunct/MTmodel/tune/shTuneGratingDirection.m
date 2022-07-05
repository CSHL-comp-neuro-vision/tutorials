% [xDirection, yResponse] = shTuneGratingDirection(pars, neuron, stageName, 
%                           nDataPoints, gratingSf, gratingTf, gratingContrast)
%
% shTuneGratingDirection computes a tuning curve of response vs. direction 
% for a full field drifting grating.
%
% Output:
% xDirection        a vector containing the directions of all the stimuli.
% yResponse         the average response of the neuron to each stimulus.
%
% Required arguments:
% pars              a parameters structure
% neuron            parameters of the neuron you want to model
% stageName         name of the stage you want to model. Choices:
%                   'v1Complex' and 'mtPattern'
%
% Optional arguments:
% nDataPoints       the number of different stimulus conditions. 
%                   DEFAULT = 7
% gratingSf         the spatial frequency of the grating, in cycles/pixel.
%                   DEFAULT = the neuron's preferred spatial frequency.
% gratingTf         the temporal frequency of the grating, in cycles/frame.
%                   DEFAULT = the neuron's preferred temporal
%                   frequency.
% gratingContrast   the contrast of the grating. DEFAULT = 1


function [xDirection, yResponse] = shTuneGratingDirection(varargin)

% The following arguments are optional and are by default 'default'.
nDataPoints = 'default';
gratingSf = 'default';
gratingTf = 'default';
gratingContrast = 'default';

% parse the varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};              end
if nargin >= 5;     gratingSf = varargin{5};                end
if nargin >= 6;     gratingTf = varargin{6};                end
if nargin >= 7;     gratingContrast = varargin{7};          end

if strcmp(stageName(1:2), 'v1');
    preferredGrating = v12sin(neuron);
else
    preferredGrating = mt2sin(neuron);
end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 7;                        end
if strcmp(gratingSf, 'default');            gratingSf = preferredGrating(2);        end
if strcmp(gratingTf, 'default');            gratingTf = preferredGrating(3);        end
if strcmp(gratingContrast, 'default');      gratingContrast = 1;                    end

% Done parsing arguments. Now on with the code.
stimSz = shGetDims(pars, stageName, [1 1 31]);
xDirection = linspace(0, 2.*pi, nDataPoints);
yResponse = zeros(1, length(xDirection));
for i = 1:length(xDirection)

    % Make the stimulus
    thisStimulus = mkSin(stimSz, xDirection(i), gratingSf, gratingTf, ...
                         gratingContrast);
    
    % get the neuron's response to this stimulus
    [pop, ind, res] = shModel(thisStimulus, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
    yResponse(i) = mean(shGetNeuron(res, ind));

    % plot the results so far
    plot(180*xDirection(1:i)/pi, yResponse(1:i), 'r-', 180*xDirection(1:i)/pi, yResponse(1:i), 'k.');
    xlabel('direction (deg)'); ylabel('response');
    axis([0 360 min([0, yResponse]), max(.00001, 1.2*max(yResponse))]);
    drawnow
end
