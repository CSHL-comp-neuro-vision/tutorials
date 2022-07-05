% [xCoherence, yResponse] = shTuneDotCoherence(pars, neuron, stageName,
%                              nDataPoints, nFrames, dotDirection, dotSpeed, 
%                              dotDensity, dotSize, )
%
% shTuneDotCoherence returns a tuning curve of response to drifting dot
% stimulus vs. motion coherence. 
%
% Output:
% xCoherence        the coherence values at which the response is measured
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
% dotDensity        The density of the dots. DEFAULT = the default density
%                   set by mkDots.
% dotRadius         The size of the dots. If dotRadius < 0, pixel dots are
%                   used. If dotSz > 0, gaussian dots with sigma = .5 *
%                   dotRadius are used. DEFAULT = -1.
% dotPlacementStyle see documentation in mkDots for explanation. DEFAULT =
%                   'exact'

function [xCoherence, yResponse] = shTuneDotCoherence(varargin)

% By default, the optional arguments are 'default'
nDataPoints = 'default';
nFrames = 'default';
dotDirection = 'default';
dotSpeed = 'default';
dotDensity = 'default';
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
if nargin >= 8;     dotDensity = varargin{8};           end
if nargin >= 9;     dotRadius = varargin{9};            end
if nargin >= 10;    dotPlacementStyle = varargin{10};   end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 7;                end
if strcmp(nFrames, 'default');              nFrames = 71;                   end
if strcmp(dotDirection, 'default');         dotDirection = neuron(1);       end
if strcmp(dotSpeed, 'default');             dotSpeed = neuron(2);           end
if strcmp(dotDensity, 'default');           dotDensity = 'default';         end
if strcmp(dotRadius, 'default');            dotRadius = -1;                 end
if strcmp(dotPlacementStyle, 'default');    dotPlacementStyle = 'exact';    end

% We're done parsing arguments. Now on with the code:
stimSz = shGetDims(pars, stageName, [1 1 nFrames]);
xCoherence = linspace(-(nDataPoints-1), 0, nDataPoints);
xCoherence = 2.^xCoherence;
yResponse = zeros(1, length(xCoherence));
for i = 1:length(xCoherence)

    % make stimulus
    s = mkDots(stimSz, dotDirection, dotSpeed, dotDensity, xCoherence(i), ...
                dotRadius, dotPlacementStyle);

    % calculate cell's response
    [pop, ind, res] = shModel(s, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
    yResponse(i) = mean(shGetNeuron(res, ind));

    % display results so far
    plot(xCoherence(1:i), yResponse(1:i), 'r-', xCoherence(1:i), yResponse(1:i), 'k.');
    xlabel('Stimulus coherence'); ylabel('response');
    axis([min(xCoherence) max(xCoherence) min([0, yResponse]), max(1.2*max(yResponse), .00001)]);
    drawnow
end