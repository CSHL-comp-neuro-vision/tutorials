% [xDirection, yResponse] = shTuneDotMaskDirection(pars, neuron, stageName, reso,
%                                    nFrames, dotSpeed, dotDensity,
%                                    dotCoherence, dotSize)
%
% Show two superimposed drifting dot patterns to a neuron: one moving in
% the preferred direction, and a second mask pattern whose direction varies
% in each stimulus condition. Measure the average response of the neuron
% for each direction in which the mask moves.
%
% xDirection is the direction of the mask's movement, which varies from 0 to 2*pi
% yResponse is the cell's average response for each stimulus condition.
%
% PARS is a parameters structure.
% NEURON is the tuning of the neuron you want to model in V1 or MT
% coordinates.
% stageName is the name of the stage you want to model
% RESO is the number of data points you want
% NFRAMES is the number of frames you want the stimulus to be. A higher
%   number will result in a smoother tuning curve, but a longer running time.
% DOTSPEED is the speed of the dots in pixels per frame. If DOTSPEED is not
%   supplied, or if it is the string 'default', the dot speed will be set to
%   the neuron's preferred speed.
% DOTDENSITY is the density of the dot pattern. If DOTDENSITY is not
%   supplied, or it is the string 'default', then DOTDENSITY will be .1
% DOTCOHERENCE is the coherence of the dot pattern. If DOTCOHERENCE s not
%   supplied, or if it is the string 'default', then DOTCOHERENCE will be 1.
% DOTSIZE is the size of the dots. If DOTSZ is not supplied, or if it is
%   the string 'default', the dots will be impulses.

function [xDirection, yResponse] = shTuneDotMaskDirection(varargin)

% PARSE THE ARGUMENTS
pars = varargin{1};
neuron = varargin{2};
stageName = varargin{3};
%
if nargin >= 4
    reso = varargin{4};
else
    reso = 7;
end
%
if nargin >= 5
    nFrames = varargin{5};
else
    nFrames = 'default';
end
if strcmp(nFrames, 'default')
    nFrames = 101;
end
%
if nargin >= 6
    dotSpeed = varargin{6};
else
    dotSpeed = 'default';
end
if strcmp(dotSpeed, 'default')
    dotSpeed = neuron(2);
end
%
if nargin >= 7
    dotDensity = varargin{7};
else
    dotDensity = 'default';
end
if strcmp(dotDensity, 'default')
    dotDensity = .1;
end
%
if nargin >= 8
    dotCoherence = varargin{8};
else
    dotCoherence = 'default';
end
if strcmp(dotCoherence, 'default')
    dotCoherence = 1;
end
%
if  nargin >= 9
    dotSize = varargin{9};
else
    dotSize = 'default';
end
if strcmp(dotSize, 'default')
    dotSize = -1;
end
% END PARSING ARGUMENTS

dims = shGetDims(pars, stageName, [1 1 nFrames]);
xDirection = linspace(neuron(1) - pi, neuron(1) + pi, reso);
xDirection = unique(sort(mod(xDirection, 2*pi)));
yResponse = zeros(1, length(xDirection));
for i = 1:length(xDirection)  
    % make the stimulus
    s1 = mkDots(dims, neuron(1), dotSpeed, dotDensity, dotCoherence, dotSize);
    s2 = mkDots(dims, xDirection(i), dotSpeed, dotDensity, dotCoherence, dotSize);
    s = s1 + s2;
    s(s>1) = 1;
    
    % get the neuron's average response to the stimulus
    [pop, ind, res] = shModel(s, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
    yResponse(i) = mean(shGetNeuron(res, ind));
    
    % plot the response so far
    plot(180*xDirection(1:i)/pi, yResponse(1:i), 'r-', 180*xDirection(1:i)/pi, yResponse(1:i), 'k.');
    xlabel('direction (deg)'); ylabel('response');
    axis([0 360 min([0, yResponse]), 1.2*max(yResponse)]);
    drawnow
end
