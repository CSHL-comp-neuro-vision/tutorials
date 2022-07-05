% s = mkWin(stimSz, windowSz, windowEdgeWidth, windowPosition);
%
% Make a cosine edged window.
%
% Required arguments:
% stimSz            The size of the reponse in [Y X T] coordinates
% windowSz          A 2-vector specifying size of the window. The first element
%                   is required, and specifies the radius of the window. The 
%                   second element is optional and specifies the duration of
%                   the window. If the duration isn't supplied, the window will
%                   last for the entire length of the stimulus.
%
% Optional arguments:
% windowEdgeWidth   A 2-vector specifying the width of the smooth cosine 
%                   edge of the window. The first element is the width of
%                   the spatial edge, the second element is the width of
%                   the temporal edge. DEFAULT = 1/4 of the window size.
%                   The window edges are added to the window, so if the
%                   radius of the window is 5 and the edge is 3, the total
%                   radius of the window will be 8.
% windowPosition    The center point of the window. DEFAULT = the center of
%                   the stimulus.

function s =  mkWin(varargin)

% the following arguments are optional and by default are 'default'
windowEdgeWidth = 'default';
windowPosition = 'default';

% Parse the varargin
stimSz = varargin{1};
windowSz = varargin{2};
if nargin >= 3;     windowEdgeWidth = varargin{3};  end
if nargin >= 4;     windowPosition = varargin{4};   end

% Assign default values where appropriate
if strcmp(windowEdgeWidth, 'default');  windowEdgeWidth = windowSz./4;          end;
if strcmp(windowPosition, 'default');   windowPosition = floor(stimSz./2)+1;    end;

% If user hasn't supplied a duration for the window, make the window last
% as long as the stimulus.
windowRadius = windowSz(1);
if length(windowSz) == 1;
    windowDuration = stimSz(3);
else
    windowDuration = windowSz(2);
end

windowSpatialEdgeWidth = windowEdgeWidth(1);
if length(windowEdgeWidth) == 1
    windowTimeEdgeWidth = windowDuration/4;
else
    windowTimeEdgeWidth = windowEdgeWidth(2);
end


% Now make the window. We'll construct it in two pieces: spatial and
% temporal.
[xSpatial, ySpatial] = meshgrid([1:stimSz(1)]-windowPosition(1), ...
                                [1:stimSz(2)]-windowPosition(2));
rSpatial = sqrt(xSpatial.^2 + ySpatial.^2);
wOne = find(rSpatial <= windowRadius);
wEdge = find(rSpatial > windowRadius & ...
             rSpatial < windowRadius + windowSpatialEdgeWidth);
wZero = find(rSpatial >= windowRadius + windowSpatialEdgeWidth);

spatialWindow = zeros(stimSz(1), stimSz(2));
spatialWindow(wOne) = 1;
spatialWindow(wEdge) = cos((pi/2)*(rSpatial(wEdge) - windowRadius)/windowSpatialEdgeWidth);
spatialWindow(wZero) = 0;
spatialWindow = repmat(spatialWindow, [1 1 stimSz(3)]);

% Now make the temporal window
xTemporal = abs([1:stimSz(3)]-windowPosition(3));
wOne = find(xTemporal <= windowDuration/2);
wEdge = find(xTemporal > windowDuration/2 & ...
             xTemporal < windowDuration/2 + windowTimeEdgeWidth);
wZero = find(xTemporal >= windowDuration/2 + windowTimeEdgeWidth);

temporalWindow = zeros(1, 1, stimSz(3));
temporalWindow(wOne) = 1;
temporalWindow(wEdge) = cos((pi/2)*(xTemporal(wEdge) - windowDuration/2)/windowTimeEdgeWidth);
temporalWindow(wZero) = 0;
temporalWindow = repmat(temporalWindow, [stimSz(1) stimSz(2), 1]);

% Take the pointwise product of the two windows to get the final result
s = spatialWindow .* temporalWindow;