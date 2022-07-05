% s = mkBar(stimSz, barDirection, barSpeed, barWidth, barEdgeWidth, 
%           pixelsBetweenBars, barStartingPosition);
%
% Make a drifting bar stimulus.
%
% Required arguments:
% stimSz                the size of the entire stimulus, in [Y X T] coordinates
% barDirection          the direction of the bar motion in radians with 0 = right
% barSpeed              the speed of the bar motion in frames/second
%
% Optional arguments:
% barWidth              the width of the center of the bar in pixels.
%                       DEFAULT = 1
% barEdgeWidth          the width of the cosine edges of the bar in pixels. An
%                       edge of width barEdgeWidth will be tacked onto both
%                       sides of the bar, so the total width will be 
%                       barWidth + 2*barEdgeWidth. DEFAULT = 2.
% pixelsBetweenBars     the number of pixels between the bars in the stimulus
%                       DEFAULT = stimSz(1)/4.
% barStartingPosition   the starting position of the first bar in pixels.
%                       The coordinate system is a line lying along the
%                       direction of the bar motion and passing through the
%                       center of the stimulus. The point 0 is the center
%                       of the stimulus. DEAFULT = 0

function driftingBar = mkBar(varargin)

% The following arguments are optional and by default are 'default'
barWidth = 'default';
barEdgeWidth = 'default';
pixelsBetweenBars = 'default';
barStartingPosition = 'default';

% parse the varargin
                        stimSz = varargin{1};
                        barDirection = varargin{2};
                        barSpeed = varargin{3};
if nargin >= 4;         barWidth = varargin{4};             end
if nargin >= 5;         barEdgeWidth = varargin{5};         end
if nargin >= 6;         pixelsBetweenBars = varargin{6};    end
if nargin >= 7;         barStartingPosition = varargin{7};  end

% Assign default values where appropriate
if strcmp(barWidth, 'default');             barWidth = 5;                       end
if strcmp(barEdgeWidth, 'default');         barEdgeWidth = 3;                   end
if strcmp(pixelsBetweenBars, 'default');    pixelsBetweenBars = stimSz(1)/4;    end
if strcmp(barStartingPosition, 'default');  barStartingPosition = 0;            end

% Calculate a few things.
gratingSf = 1/pixelsBetweenBars;
gratingTf = gratingSf * barSpeed;


% The plan of action: make a drifting grating with the right velocity and
% orientation. The period of the sin wave will be the same as the
% distance between the bars that the user wants. Then we make bars out of
% the peaks of the sin wave.

% Make the grating
gratingPhase = barStartingPosition * gratingSf;
driftingBar = 2*mkSin(stimSz, barDirection, gratingSf, gratingTf, 1, gratingPhase) - 1;

% Find the thresholds
barInnerThreshold = cos(2*pi*gratingSf*barWidth/2);
barOuterThreshold = cos(2*pi*gratingSf*(barWidth/2 + barEdgeWidth));

% There are three regions: where the stimulus should be one (the centers of
% the bars), where it should be zero (outside the bars), and where it
% should be somehwere in between (edges of the bars). Find them
wOne = find(driftingBar >= barInnerThreshold);
wEdge = find(driftingBar < barInnerThreshold & ...
             driftingBar > barOuterThreshold);
wZero = find(driftingBar <= barOuterThreshold);

% Set the regions to the appropriate level
driftingBar(wOne) = 1;
driftingBar(wZero) = 0;

% keyboard
driftingBar(wEdge) = acos(driftingBar(wEdge));      % now it ranges from 0 to 2*pi
driftingBar(wEdge) = driftingBar(wEdge)/(2*pi*gratingSf);       % now it ranges from 0 to 1 spatial period
driftingBar(wEdge) = (pi/2)*(driftingBar(wEdge) - barWidth/2)/(barEdgeWidth);
driftingBar(wEdge) = cos(driftingBar(wEdge));
