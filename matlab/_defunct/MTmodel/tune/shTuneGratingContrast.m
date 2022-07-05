% [xContrast, yResponse] = shTuneGratingContrast(pars, neuron, stageName, 
%                           nDataPoints, gratingDirection, gratingSf, gratingTf)
%
% shTuneGratingContrast computes a tuning curve of response vs. contrast for a
% full field drifting grating.
%
% Output:
% xContrast         a vector containing the contrasts of all the stimuli.
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
% gratingDirection  the direction of the grating, in radians, 0 =
%                   rightward. DEFAULT = the neuron'thisStimulus preferred direction.
% gratingSf         the spatial frequency of the grating, in cycles/pixel.
%                   DEFAULT = the neuron'thisStimulus preferred spatial frequency.
% gratingTf         the temporal frequency of the grating, in cycles/frame.
%                   DEFAULT = the neuron'thisStimulus preferred temporal frequency.

function [xContrast, yResponse] = shtunedir(varargin)

% The following arguments are optional and are by default 'default'.
nDataPoints = 'default';
gratingDirection = 'default';
gratingSf = 'default';
gratingTf = 'default';

% parse the varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};              end
if nargin >= 5;     gratingDirection = varargin{5};         end
if nargin >= 6;     gratingSf = varargin{6};                end
if nargin >= 7;     gratingTf = varargin{7};                end

if strcmp(stageName(1:2), 'v1');
    preferredGrating = v12sin(neuron);
else
    preferredGrating = mt2sin(neuron);
end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 7;                        end
if strcmp(gratingDirection, 'default');     gratingDirection = preferredGrating(1); end
if strcmp(gratingSf, 'default');            gratingSf = preferredGrating(2);        end
if strcmp(gratingTf, 'default');            gratingTf = preferredGrating(3);        end

% Done parsing arguments. Now on with the code.
stimSz = shGetDims(pars, stageName, [1 1 31]);
xContrast = logspace(-2, 0, nDataPoints);
yResponse = zeros(size(xContrast));
for i = 1:length(xContrast)

    % generate the stimulus    
    thisStimulus = mkSin(stimSz, gratingDirection, gratingSf, gratingTf, ...
                         xContrast(i));
    
    % get the neuron's response to this stimulus
    [pop, ind, res] = shModel(thisStimulus, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
    yResponse(i) = mean(shGetNeuron(res, ind));

    % plot the response so far
    semilogx(xContrast(1:i), yResponse(1:i), 'r-', xContrast(1:i), yResponse(1:i), 'k.');
    xlabel('quote unquote contrast'); ylabel('response');
    axis([min(xContrast) max(xContrast) 0 max(.000001, 1.2.*max(yResponse))]);
    drawnow;
end
