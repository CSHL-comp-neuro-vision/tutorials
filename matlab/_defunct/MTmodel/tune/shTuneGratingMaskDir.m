% [xDirection, yResponse] = shTuneGratingMaskDir(pars, neuron, stageName, 
%                           nDataPoints, gratingDirection, 
%                           gratingSf, gratingTf, stimulusContrast)
%
% shTuneGratingMaskDir computes a tuning curve of response vs. direction of
% mask for a stimuli composed of a grating + a superimposed mask.
%
% Output:
% xMaskDirection    a vector containing the directions of the mask in each
%                   stimulus condition. 
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
%                   rightward. DEFAULT = the neuron'thisStimulus preferred
%                   direction.
% gratingSf         the spatial frequency of the grating, in cycles/pixel.
%                   DEFAULT = the neuron's preferred spatial frequency.
% gratingTf         the temporal frequency of the grating, in cycles/frame.
%                   DEFAULT = the neuron's preferred temporal
%                   frequency.
% stimulusContrast  the contrast of the stimulus. DEFAULT = 1


function [xMaskDirection, yResponse] = shTuneGratingMaskDir(varargin)

% The following arguments are optional and are by default 'default'.
nDataPoints = 'default';
gratingDirection = 'default';
gratingSf = 'default';
gratingTf = 'default';
stimulusContrast = 'default';

% parse the varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};              end
if nargin >= 5;     gratingDirection = varargin{5};         end
if nargin >= 6;     gratingSf = varargin{6};                end
if nargin >= 7;     gratingTf = varargin{7};                end
if nargin >= 8;     stimulusContrast = varargin{8};         end

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
if strcmp(stimulusContrast, 'default');     stimulusContrast = 1;                   end


% Done parsing arguments. Now get on with it.
stimSz = shGetDims(pars, stageName, [1 1 31]);
xMaskDirection = linspace(0, 2*pi, nDataPoints);
yResponse = zeros(1, length(xMaskDirection));
for i = 1:length(xMaskDirection)

    % generate the stimulus
    thisGrating = mkSin(stimSz, gratingDirection, gratingSf, gratingTf, ...
                        stimulusContrast/2);
    thisMaskGrating = mkSin(stimSz, xMaskDirection(i), gratingSf, ...
                            gratingTf, stimulusContrast/2);
    thisStimulus = thisGrating + thisMaskGrating;

    % Calculate the neuron's response to the stimulus.
    [pop, ind, res] = shModel(thisStimulus, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
    yResponse(i) = mean2(shGetNeuron(res, ind, 1, 1));
    
    % plot the results so far
    plot(180*xMaskDirection(1:i)/pi, yResponse(1:i), 'r-', 180*xMaskDirection(1:i)/pi, yResponse(1:i), 'k.');
    xlabel('direction (deg)'); ylabel('response');
    axis([0 360 min([0, yResponse]), 1.2*max(yResponse)]);
    drawnow
end
