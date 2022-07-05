% [xDirection, yResponse] = shTuneDotDirection(pars, neuron, stageName,
%                              nDataPoints, nFrames, dotSpeed, dotDensity, 
%                              dotCoherence, dotRadius, dotPlacementStyle)
%
% shTuneDotDirection returns a tuning curve of response to a drifting dot
% stimulus vs. direction of motion.
%
% Output:
% xDirection        the direction of motion for each stimulus
% yResponse         the response for each entry in xCoherence
%
% Required arguments:
% pars              parameters structure for the model.
% neuron            tuning of the neuron you want to model
% stageName         name of the model stage you want. Options 'v1Complex',
%                   'mtPattern'.
%
% Optional arguments:
% nFrames           The number of frames in each stimulus. The higher
%                   nFrames is, the smoother and more regular the tuning
%                   curve will be, but the longer the code will take to
%                   run. DEFAULT = 71;
% dotSpeed          The speed of the dot motion in pixels/frame. DEFAULT =
%                   the neuron's preferred speed.
% dotDensity        The density of the dots. DEFAULT = the default density
%                   set by mkDots.
% dotCoherence      The coherence of dot motion. DEFAULT = 1.
% dotRadius         The size of the dots. If dotRadius < 0, pixel dots are
%                   used. If dotSz > 0, gaussian dots with sigma = .5 *
%                   dotRadius are used. DEFAULT = -1.
% dotPlacementStyle see documentation in mkDots for explanation. DEFAULT =
%                   'exact'

function [xDirection, yResponse] = shTuneDotDirection(varargin);

% By default, the optional arguments are 'default'
nDataPoints = 'default';
nFrames = 'default';
dotSpeed = 'default';
dotDensity = 'default';
dotCoherence = 'default';
dotRadius = 'default';
dotPlacementStyle = 'default';

% Parse arguments out of varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};          end
if nargin >= 5;     nFrames = varargin{5};              end
if nargin >= 6;     dotSpeed = varargin{6};             end
if nargin >= 7;     dotDensity = varargin{7};           end
if nargin >= 8;     dotCoherence = varargin{8};         end
if nargin >= 9;     dotRadius = varargin{9};            end
if nargin >= 10;    dotPlacementStyle = varargin{10};   end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 9;                end
if strcmp(nFrames, 'default');              nFrames = 71;                   end
if strcmp(dotSpeed, 'default');             dotSpeed = neuron(2);           end
if strcmp(dotDensity, 'default');           dotDensity = 'default';         end
if strcmp(dotCoherence, 'default');         dotCoherence = 1;               end
if strcmp(dotRadius, 'default');            dotRadius = -1;                 end
if strcmp(dotPlacementStyle, 'default');    dotPlacementStyle = 'exact';    end

% We're done parsing arguments. Now on with the code.
stimSz = shGetDims(pars, stageName, [1 1 nFrames]);
xDirection = linspace(neuron(1) - pi, neuron(1) + pi, nDataPoints+1);
xDirection = unique(sort(mod(xDirection, 2*pi)));
yResponse = zeros(1, length(xDirection));
for i = 1:length(xDirection)

    % make the stimulus
    s = mkDots(stimSz, xDirection(i), dotSpeed, dotDensity, dotCoherence, ...
        dotRadius, dotPlacementStyle);

    % calculate the neuron's response
    [pop, ind, res] = shModel(s, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
    yResponse(i) = mean2(shGetNeuron(res, ind, 1, 1));

    % plot the results so far
    plot(180*xDirection(1:i)/pi, yResponse(1:i), 'r-', 180*xDirection(1:i)/pi, yResponse(1:i), 'k.');
    xlabel('direction (deg)'); ylabel('response');
    axis([0 360 min([0, yResponse]), max(.00005, 1.2*max(yResponse))]);
    drawnow
end