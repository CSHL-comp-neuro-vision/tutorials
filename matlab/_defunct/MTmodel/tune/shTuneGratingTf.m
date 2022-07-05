% [xTf, yResponse] = shTuneGratingTf(pars, neuron, stageName, nDataPoints,
%                                    gratingDirection, gratingSf, gratingContrast)
%
% shTuneGratingTf computes a tuning curve of response vs. temporal frequency 
% for a full field drifting grating.
%
% Output:
% xTf               a vector containing the temporal frequencies of all the stimuli.
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
% gratingDirection  the direction of the nonmask grating, in radians, 0 =
%                   rightward. DEFAULT = the neuron's preferred
%                   direction.
% gratingSf         the spatial frequency of the grating, in cycles/pixel.
%                   DEFAULT = the neuron's preferred spatial
%                   frequency.
% gratingContrast   the contrast of the grating. DEFAULT = 1


function [xTf, yResponse] = shTuneGratingTf(varargin)

% The following arguments are optional and are by default 'default'.
nDataPoints = 'default';
gratingDirection = 'default';
gratingSf = 'default';
gratingContrast = 'default';

% parse the varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};              end
if nargin >= 5;     gratingDirection = varargin{5};         end
if nargin >= 6;     gratingSf = varargin{6};                end
if nargin >= 7;     gratingContrast = varargin{7};          end

if strcmp(stageName(1:2), 'v1');
    preferredGrating = v12sin(neuron);
else
    preferredGrating = mt2sin(neuron);
end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 7;                        end
if strcmp(gratingDirection, 'default');     gratingDirection = preferredGrating(1); end
if strcmp(gratingSf, 'default');            gratingSf = preferredGrating(3);        end
if strcmp(gratingContrast, 'default');      gratingContrast = 1;                    end

% done parsing arguments. Now get on with it.
stimSz = shGetDims(pars, stageName, [1 1 31]);
xTf = logspace(-2, -.3010, nDataPoints);
yResponse = zeros(size(xTf));
for i = 1:length(xTf)

    % generate this stimulus
    thisStimulus = mkSin(stimSz, gratingDirection, gratingSf, xTf(i), gratingContrast);
    [pop, ind, res] = shModel(thisStimulus, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        tmp = sqrt(res.^2);
    end
    yResponse(i) = mean(shGetNeuron(res, ind));

    % plot the results thus far
    semilogx(xTf(1:i), yResponse(1:i), 'r-', xTf(1:i), yResponse(1:i), 'k.');
    xlabel('sf (c/px)'); ylabel('response');
    axis([min(xTf) max(xTf) 0 max(.000005, 1.2*max(yResponse))]);
    drawnow
end
