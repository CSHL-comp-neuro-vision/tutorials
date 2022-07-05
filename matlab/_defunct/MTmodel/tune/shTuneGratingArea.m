% [xRadius, yResponse] = shTuneGratingArea(pars, neuron, stageName, nDataPoints,
%                               minRadius, maxRadius, gratingDirection, 
%                               gratingSf, gratingTf, gratingContrast)
%
% shTuneGratingArea computes a tuning curve of response vs. radius for a
% windowed drifting grating.
%
% Output:
% xRadius           a vector of the radius of the window in each stimulus.
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
% minRadius         the smallest radius to test, in pixels. DEFAULT = 1.
% maxRadius         the largest radius to test, in pixels. DEFAULT = 15.
% gratingDirection  the direction of the grating, in radians, 0 =
%                   rightward. DEFAULT = the neuron's preferred direction.
% gratingSf         the spatial frequency of the grating, in cycles/pixel.
%                   DEFAULT = the neuron's preferred spatial frequency.
% gratingTf         the temporal frequency of the grating, in cycles/frame.
%                   DEFAULT = the neuron's preferred temporal frequency.
% gratingContrast   the contrast of the grating. DEFAULT = 1


function [xRadius, yResponse] = shTuneGratingArea(varargin)

% The following arguments are optional and are by default 'default'.
nDataPoints = 'default';
minRadius = 'default';
maxRadius = 'default';
gratingDirection = 'default';
gratingSf = 'default';
gratingTf = 'default';
gratingContrast = 'default';

% parse the varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};              end
if nargin >= 5;     minRadius = varargin{5};                end
if nargin >= 6;     maxRadius = varargin{6};                end
if nargin >= 7;     gratingDirection = varargin{7};         end
if nargin >= 8;     gratingSf = varargin{8};                end
if nargin >= 9;     gratingTf = varargin{9};                end
if nargin >= 10;    gratingContrast = varargin{10};         end

if strcmp(stageName(1:2), 'v1');
    preferredGrating = v12sin(neuron);
else
    preferredGrating = mt2sin(neuron);
end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 7;                        end
if strcmp(minRadius, 'default');            minRadius = 1;                          end
if strcmp(maxRadius, 'default');            maxRadius = 15;                         end
if strcmp(gratingDirection, 'default');     gratingDirection = preferredGrating(1); end
if strcmp(gratingSf, 'default');            gratingSf = preferredGrating(2);        end
if strcmp(gratingTf, 'default');            gratingTf = preferredGrating(3);        end
if strcmp(gratingContrast, 'default');      gratingContrast = 1;                    end

% Done parsing arguments. Get on with it.
stimSz = shGetDims(pars, stageName, [1 1 31]);
xRadius = linspace(minRadius, maxRadius, nDataPoints);
yResponse = zeros(length(xRadius), 1);
for i = 1:length(xRadius)

    % Make the stimulus
    stimulusGrating = mkSin(stimSz, gratingDirection, gratingSf, gratingTf, ...
                            gratingContrast);
    windowEdgeWidth = 2;
    stimulusWindow = mkWin(stimSz, xRadius(i), windowEdgeWidth);
    thisStimulus = stimulusGrating .* stimulusWindow;

    % Calculate the neuron's response to this stimulus
    [pop, ind, res] = shModel(thisStimulus, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end

    % Plot the results so far.
    yResponse(i) = mean(shGetNeuron(res, ind));
    plot(xRadius(1:i), yResponse(1:i), 'r-', xRadius(1:i), yResponse(1:i), 'k.');
    xlabel('radius (px)'); ylabel('response');
    axis([0, xRadius(end), 0, 1.2.*max(yResponse)]);
    drawnow
end
