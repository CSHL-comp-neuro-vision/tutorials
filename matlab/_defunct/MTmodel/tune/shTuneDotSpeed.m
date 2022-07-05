% [xSpeed, yResponse] = shTuneDotDensity(pars, neuron, stageName,
%                              nDataPoints, nFrames, dotDirection, dotDensity, 
%                              dotCoherence, dotRadius, dotPlacementStyle)
%
% shTuneDotSpeed returns a tuning curve of response to a drifting dot
% stimulus vs. speed of motion.
%
% Output:
% xSpeed            the speed of motion for each stimulus
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
% dotDirection      The direction of the dot motion in pixels/frame. DEFAULT =
%                   the neuron's preferred direction.
% dotDensity        The density of the dots. DEFAULT = the default density
%                   set by mkDots.
% dotCoherence      The coherence of dot motion. DEFAULT = 1.
% dotRadius         The size of the dots. If dotRadius < 0, pixel dots are
%                   used. If dotSz > 0, gaussian dots with sigma = .5 *
%                   dotRadius are used. DEFAULT = -1.
% dotPlacementStyle see documentation in mkDots for explanation. DEFAULT =
%                   'exact'


function [xSpeed, yResponse] = shTuneDotSpeed(varargin)

% By default, the optional arguments are 'default'
nDataPoints = 'default';
nFrames = 'default';
dotDirection = 'default';
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
if nargin >= 6;     dotDirection = varargin{6};         end
if nargin >= 7;     dotDensity = varargin{7};           end
if nargin >= 8;     dotCoherence = varargin{8};         end
if nargin >= 9;     dotRadius = varargin{9};            end
if nargin >= 10;    dotPlacementStyle = varargin{10};   end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 9;                end
if strcmp(nFrames, 'default');              nFrames = 71;                   end
if strcmp(dotDirection, 'default');         dotDirection = neuron(1);       end
if strcmp(dotDensity, 'default');           dotDensity = 'default';         end
if strcmp(dotCoherence, 'default');         dotCoherence = 1;               end
if strcmp(dotRadius, 'default');            dotRadius = 3;                  end
if strcmp(dotPlacementStyle, 'default');    dotPlacementStyle = 'exact';    end

% Done parsing arguments. Now get to work.
stimSz = shGetDims(pars, stageName, [25 25 nFrames]);
prefSpeed = neuron(2);
minSpeed = log2(.25.*prefSpeed);
maxSpeed = log2(4.*prefSpeed);

xSpeed = linspace(minSpeed, maxSpeed, nDataPoints);
xSpeed = 2.^xSpeed;
yResponse = zeros(1, length(xSpeed));

for i = 1:length(xSpeed)
    % make the stimulus
    s = mkDots(stimSz, dotDirection, xSpeed(i), dotDensity, dotCoherence, ...
                dotRadius, dotPlacementStyle);

%     % apply a lowpass filter
%     f = fftshift(fftn(ifftshift(s1)));
%     centerPoint = floor(stimSz./2) + 1;
%     fWindow = mkWin(stimSz, [7 31], [3 11], centerPoint);
%     f = f.*fWindow;
%     s = fftshift(ifftn(ifftshift(f)));
%     s = real(s);

    % get the neuron's response to the stimulus
    [pop, ind, res] = shModel(s, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
    yResponse(i) = mean2(shGetNeuron(res, ind));

    % plot the results thus far
    semilogx(xSpeed(1:i), yResponse(1:i), 'r-', xSpeed(1:i), yResponse(1:i), 'k.');
    xlabel('speed (px/frame)'); ylabel('response');
    axis([min(xSpeed) max(xSpeed) min([0, yResponse]), max(.00001, 1.2*max(yResponse))]);
    drawnow
end
