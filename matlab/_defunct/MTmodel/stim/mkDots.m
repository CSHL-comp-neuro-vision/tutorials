% function s = mkDots(stimSz, dotDirection, dotSpeed, dotDensity, dotCoherence,
%                   dotRadius, densityStyle, sampleFactor, interpMethod)
%
% MKDOTS makes a drifting dot stimulus.
%
% Required arguments:
% stimSz            The dimensions of the stimulus, in [Y X T] coordinates;
% dotDirection      The direction of movement, in radians, with 0 = rightward.
%                   If dotDirection is just one number, the dots will move
%                   in that direction at all times. If dotDirection is a
%                   vector of the same length as the number of frames in
%                   the stimulus, the direction of movement in each frame
%                   will be specified by the corresponding element of
%                   dotDirection.
% dotSpeed          The speed of movement, in frames/second.
%                   If dotSpeed is just one number, the dots will move
%                   with that speed at all times. If dotSpeed is a
%                   vector of the same length as the number of frames in
%                   the stimulus, the speed of movement in each frame
%                   will be specified by the corresponding element of
%                   dotSpeed.
%
% Optional arguments:
% dotDensity        The density of the dots, which can be from 0 to 1. DEFAULT = .1
% dotCoherence      The coherence of dot motion. DEFAULT = 1.
% dotRadius         The radius of the dots. If dotRadius < 0, the dots will
%                   be single pixels. If dotRadius > 0, the dots will be
%                   Gaussian blobs with sigma = .5 * dotRadius. DEFAULT = -1                    
% dotPlacementStyle The number of dots is calculated by multiplying the
%                   dotDensity by the size of the image window. If the dots
%                   are placed randomly, some dots will overlap, so the
%                   actual dot density will be lower than dotDensity. If,
%                   however, the dots are placed exactly so that they never
%                   overlap, all the dots will move in register. Both of
%                   these are problems that only occur when you use pixel
%                   dots rather than Gaussian dots. Your choices for
%                   dotPlacementStyle are 'random' (for the first problem)
%                   or 'exact' (for the second problem). DEFAULT = 'random'
% sampleFactor      Only important if you are using Gaussian dots. The dots
%                   are made by calculating one large Gaussian dot, and
%                   interpolating from that dot to make the smaller dots.
%                   This parameter specifies how much larger the large dot
%                   is than the smaller dots. DEFAULT = 10.
% interpMethod      Only important if you are using Gaussian dots.
%                   Specifies the interpolation method used in creating the
%                   smaller dots. Choices: 'linear', 'cubic', 'nearest',
%                   'v4'. DEFAULT = 'linear'

function s = mkDots(varargin)

% By default, these variables will have their default values.
dotDensity = 'default';
dotCoherence = 'default';
dotRadius = 'default';
densityStyle = 'default';
sampleFactor = 'default';
interpMethod = 'default';

% Make sure the user has supplied the first three arguments.
if nargin < 3
    error('You must supply at least three arguments: stimulus size, direction, and speed');
end

% Parse arguments out of varargin
stimSz = varargin{1};
dotDirection = varargin{2};
dotSpeed = varargin{3};
if nargin >= 4; dotDensity = varargin{4};       end
if nargin >= 5; dotCoherence = varargin{5};     end
if nargin >= 6; dotRadius = varargin{6};        end
if nargin >= 7; densityStyle = varargin{7};     end
if nargin >= 8; sampleFactor = varargin{8};     end
if nargin >= 9; interpMethod = varargin{9};     end

% set default values
if strcmp(dotCoherence, 'default'); dotCoherence = 1;           end
if strcmp(dotRadius, 'default');    dotRadius = -1;             end
if strcmp(dotDensity, 'default');   
    if dotRadius < 0
        dotDensity = .1;       
    else
        dotDensity = .3./(pi*dotRadius^2);
    end
end
if strcmp(densityStyle, 'default'); densityStyle = 'rough';     end
if strcmp(sampleFactor, 'default'); sampleFactor = 10;          end
if strcmp(interpMethod, 'default'); interpMethod = '*linear';   end

% resize dotDirection and dotSpeed if necessary
if length(dotDirection) == 1
    dotDirection = repmat(dotDirection, stimSz(3), 1);
end
if length(dotSpeed) == 1
    dotSpeed = repmat(dotSpeed, stimSz(3), 1);
end

if length(dotDirection) ~= stimSz(3)
    error('If dotDirection is a vector, it must have the same number of entries as there are frames in the stimulus.');
end
if length(dotSpeed) ~= stimSz(3)
    error('If dotSpeed is a vector, it must have the same number of entries as there are frames in the stimulus.');
end

%%%%%%%%%%%%%%%%% NOW, ON WITH THE CODE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%


% make the large dot for sampling
xLargeDot = linspace(-(3/2)*dotRadius, (3/2)*dotRadius, ceil(3*dotRadius*sampleFactor));
sigmaLargeDot = dotRadius/2;
largeDot = exp(-xLargeDot.^2./(2.*sigmaLargeDot^2));
largeDot = largeDot'*largeDot;
[xLargeDot, yLargeDot] = meshgrid(xLargeDot, xLargeDot);

% There is a buffer area around the image so that we don't have to worry
% about getting wraparounds exactly right.  This buffer is twice the size
% of a dot diameter.
if dotRadius > 0
    bufferSize = 2*ceil(dotRadius*3);
else
    bufferSize = 5;
end

% the 'frame' is the field across which the dots drift. The final ouput of
% this function will consist of 'snapshots' of this frame without buffer.
% We store the size of the frame and a coordinate system for it here:
frameSzX = stimSz(2) + 2.*bufferSize;
frameSzY = stimSz(1) + 2.*bufferSize;
[xFrame, yFrame] = meshgrid([1:frameSzX], [1:frameSzY]);
yFrame = flipud(yFrame);


% nDots is the number of coherently moving dots in the stimulus.
% nDotsNonCoherent is the number of noncoherently moving dots.
nDots = round(dotCoherence.*dotDensity.*prod(size(xFrame)));
nDotsNonCoherent = round((1-dotCoherence).*dotDensity.*prod(size(xFrame)));

% Set the initial dot positions.
% dotPositions is a matrix of positions of the coherently moving dots in
%   [y, x] coordinates; each row in dotPositions stores the position of one
%   dot.
densityStyle = 'exact';
if strcmp(densityStyle, 'exact')
    z = zeros(prod(size(xFrame)), 1);
    z(1:nDots) = 1;
    ord = rand(size(z));
    [blah, ord] = sort(ord);
    z = z(ord);
    dotPositions = [yFrame(z==1), xFrame(z==1)];
else
    dotPositions = rand(nDots, 2) * [frameSzY 0; 0 frameSzX];
end

% s will store the output. After looping over each frame, we will trim away
% the buffer from s to obtain the final result.
s = zeros(frameSzY, frameSzX, stimSz(3));
toInterpolate = [-floor((3/2)*dotRadius):floor((3/2)*dotRadius)];
dSz = floor((3/2)*dotRadius);
for t = 1:stimSz(3)

    % move the positions of all the dots
    dotVelocity = [sin(dotDirection(t)), cos(dotDirection(t))];
    dotVelocity = dotVelocity*dotSpeed(t);
    dotPositions = dotPositions + repmat(dotVelocity, size(dotPositions, 1), 1);

    % wrap around for all dots that have gone past the image borders
    w = find(dotPositions(:,1) > frameSzY + .5);
    dotPositions(w,1) = dotPositions(w,1) - frameSzY;

    w = find(dotPositions(:,1) < .5);
    dotPositions(w,1) = dotPositions(w,1) + frameSzY;

    w = find(dotPositions(:,2) > frameSzX + .5);
    dotPositions(w,2) = dotPositions(w,2) - frameSzX;

    w = find(dotPositions(:,2) < .5);
    dotPositions(w,2) = dotPositions(w,2) + frameSzY;

    % add noncoherent dots and make a vector of dot positions for this
    % frame only.
    dotPositionsNonCoherent = rand(nDotsNonCoherent, 2) * [frameSzY-1 0; 0 frameSzX-1] + .5;

    % create a temporary matrix of positions for dots to be shown in this
    % frame.
    tmpDotPositions = [dotPositions; dotPositionsNonCoherent];

    % prepare a matrix of zeros for this frame
    thisFrame = zeros(size(xFrame));
    if dotRadius > 0
        % in each frame, don't show dots near the edges of the frame. That's
        % why we have a buffer. The reason we don't show them is that we don't
        % want to deal with edge handling.
        w1 = find(tmpDotPositions(:,1) > frameSzY - bufferSize + (3/2)*dotRadius);
        w2 = find(tmpDotPositions(:,1) < bufferSize - (3/2)*dotRadius);
        w3 = find(tmpDotPositions(:,2) > frameSzX - bufferSize + (3/2)*dotRadius);
        w4 = find(tmpDotPositions(:,2) < bufferSize - (3/2)*dotRadius);
        w = [w1; w2; w3; w4];
        tmpDotPositions(w, :) = [];

        % add the dots to thisFrame
        for p = 1:size(tmpDotPositions, 1)

            % find the center point of the current dot, in thisFrame
            % coordinates. This is where the dot will be placed.
            cpY = round(tmpDotPositions(p, 1));
            cpX = round(tmpDotPositions(p, 2));

            xToInterpolate = toInterpolate + (round(tmpDotPositions(p,2)) - tmpDotPositions(p,2));
            yToInterpolate = toInterpolate + (round(tmpDotPositions(p,1)) - tmpDotPositions(p,1));
            [xToInterpolate, yToInterpolate] = meshgrid(xToInterpolate, yToInterpolate);
            thisSmallDot = interp2(xLargeDot, yLargeDot, largeDot, ...
                xToInterpolate, yToInterpolate, interpMethod);

            % now add this small dot to the frame.
            thisFrame(cpY-dSz:cpY+dSz, cpX-dSz:cpX+dSz) = ...
                thisFrame(cpY-dSz:cpY+dSz, cpX-dSz:cpX+dSz) + ...
                thisSmallDot;

        end
    else
        tmpDotPositions(tmpDotPositions(:,1) > frameSzY, :) = [];
        tmpDotPositions(tmpDotPositions(:,1) < 1, :) = [];
        tmpDotPositions(tmpDotPositions(:,2) > frameSzX, :) = [];
        tmpDotPositions(tmpDotPositions(:,2) < 1, :) = [];
        tmpDotPositions = round(tmpDotPositions);

        w = sub2ind(size(thisFrame), tmpDotPositions(:,1), tmpDotPositions(:,2));

        thisFrame(w) = 1;
    end
    % Add this frame to the final output
    s(:,:,t) = flipud(thisFrame);
end
% Now trim away the buff
s = s(bufferSize+1:end-bufferSize, bufferSize+1:end-bufferSize, :);
s(s>1) = 1;