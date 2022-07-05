% [xDensity, yResponse] = shTuneDotDensity(pars, neuron, stageName,
%                              nDataPoints, nFrames, dotDirection, dotSpeed, 
%                              dotCoherence, dotSize, dotPlacementStyle)
%
% shTuneDotDensity returns a tuning curve of response to a drifting dot
% stimulus vs. dot density.
%
% Output:
% xDensity          the density values at which the response is measured
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
% dotDirection      The direction of the dot motion in radians, with 0 = 
%                   rightward motion. DEFAULT = the neuron's preferred
%                   direction of motion.
% dotSpeed          The speed of the dot motion in pixels/frame. DEFAULT =
%                   the neuron's preferred speed.
% dotCoherence      The coherence of dot motion. DEFAULT = 1.
% dotRadius         The size of the dots. If dotRadius < 0, pixel dots are
%                   used. If dotSz > 0, gaussian dots with sigma = .5 *
%                   dotRadius are used. DEFAULT = -1.
% dotPlacementStyle see documentation in mkDots for explanation. DEFAULT =
%                   'exact'


function [xDensity, yResponse] = shTuneDotDensity(varargin)

% By default, the optional arguments are 'default'
nDataPoints = 'default';
nFrames = 'default';
dotDirection = 'default';
dotSpeed = 'default';
dotCoherence = 'default';
dotRadius = 'default';
dotPlacementStyle = 'default';

% Parse arguments out of varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};          end
if nargin >= 5;     nFrames = varargin{5};              end
if nargin >= 6;     dotDirection = varargin{6};         end
if nargin >= 7;     dotSpeed = varargin{7};             end
if nargin >= 8;     dotCoherence = varargin{8};         end
if nargin >= 9;     dotRadius = varargin{9};            end
if nargin >= 10;    dotPlacementStyle = varargin{10};   end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 7;                end
if strcmp(nFrames, 'default');              nFrames = 71;                   end
if strcmp(dotDirection, 'default');         dotDirection = neuron(1);       end
if strcmp(dotSpeed, 'default');             dotSpeed = neuron(2);           end
if strcmp(dotCoherence, 'default');         dotCoherence = 1;               end
if strcmp(dotRadius, 'default');            dotRadius = -1;                 end
if strcmp(dotPlacementStyle, 'default');    dotPlacementStyle = 'exact';    end

% Parsing arguments is done. Now get to work.
stimSz = shGetDims(pars, stageName, [1 1 nFrames]);
xDensity = linspace(-7, -1, nDataPoints);
xDensity = 2.^xDensity;
yResponse = zeros(1, length(xDensity));
for i = 1:length(xDensity)
    % generate the stimulus
    s = mkDots(stimSz, dotDirection, dotSpeed, xDensity(i), dotCoherence, ...
        dotRadius, dotPlacementStyle);

    % find the neuron's response to the stimulus
    [pop, ind, res] = shModel(s, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin');
        tmp = sqrt(tmp.^2);
    end
    yResponse(i) = mean(shGetNeuron(res, ind));

    % plot the results so far
    plot(xDensity(1:i), yResponse(1:i), 'r-', xDensity(1:i), yResponse(1:i), 'k.');
    xlabel('Dot density'); ylabel('response');
    axis([min(xDensity) max(xDensity) min([0, yResponse]), max(1.2*max(yResponse), .00001)]);
    drawnow
end