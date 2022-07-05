% p = drawPlane(planeNormal, planePosition, planeSize, patchArgs)
%
% Draw a plane in the current axes. drawPlane works by drawing a square in 
% the YZ plane and then rotating it to the direction specified by the user. 
%
% Required arguments:
% planeNormal           The vector normal to the plane in [Y X T] coordinates
% 
% Optional arguments:
% planePosition         The 'center point' of the plane. While a plane has
%                       no center point, drawPlane actually draws a finite
%                       rectangle that lies on a plane. planePosition
%                       defines the center of that rectangle.
%                       DEFAULT: the center point of the current axes as
%                       found from get(gca, 'xpos'), etc.
% planeSize             The size of the rectangle to be drawn. If planeSize
%                       a 1-vector, then the rectangle will be a square. If
%                       planeSize is a 2-vector, then the first element
%                       will specify the size along the Y-axis before
%                       rotation of the plane, while the second element
%                       will specify the size along the Z-axis before
%                       rotation of the plane. DEFAULT: the size of the
%                       current axes given by get(gca, 'yLim') etc.
% patchArgs             as few or asmany arguments as you like that will be
%                       passed on to the patch command that ultimately
%                       draws the plane.
%
% Output:
% p                     the handle for the patch objects

function p = drawPlane(varargin)

planePosition = 'default';
planeSize = 'default';
patchArgs = 'default';

planeNormal = varargin{1};
if nargin >= 2;         planePosition = varargin{2};            end;
if nargin >= 3;         planeSize = varargin{3};                end;
if nargin >= 4;         patchArgs = varargin(4:end);            end;

xAxesLim = get(gca, 'xlim');
yAxesLim = get(gca, 'ylim');
zAxesLim = get(gca, 'zlim');
if strcmp(planePosition, 'default');
    clear planePosition
    planePosition(2) = xAxesLim(1) + (xAxesLim(2) - xAxesLim(1))/2;
    planePosition(1) = yAxesLim(1) + (yAxesLim(2) - yAxesLim(1))/2;
    planePosition(3) = zAxesLim(1) + (zAxesLim(2) - zAxesLim(1))/2;
end

if strcmp(planeSize, 'default');
    clear planeSize
    planeSize(1) = (yAxesLim(2) - yAxesLim(1));
    planeSize(2) = (zAxesLim(2) - zAxesLim(1));
end
if length(planeSize) == 1
    planeSize = [planeSize, planeSize];
end


if strcmp(patchArgs, 'default');        patchArgs = {};         end;

% Get the initial coordinates of the corners, and store them in a matrix
x = [0, 0, 0, 0];
y = [-planeSize(1)/2, -planeSize(1)/2, planeSize(1)/2, planeSize(1)/2];
z = [-planeSize(2)/2, planeSize(2)/2, planeSize(2)/2, -planeSize(2)/2];
xyz = [x; y; z];

% Rotate the coordinates of the corners.
planeNormalSphericalCoordinates = rec2sphere(planeNormal);
azAngle = planeNormalSphericalCoordinates(1); 
elAngle = planeNormalSphericalCoordinates(2);

azRotationMatrix = [cos(azAngle) -sin(azAngle) 0; sin(azAngle) cos(azAngle) 0; 0 0 1];
elRotationMatrix = [cos(elAngle) 0 -sin(elAngle); 0 1 0; sin(elAngle) 0 cos(elAngle)]; 

xyz = (azRotationMatrix*elRotationMatrix*xyz)';

% parse out the corner coordinates and shift them to the position specified
% by planePosition
x = xyz(:,1) + planePosition(2);
y = xyz(:,2) + planePosition(1);
z = xyz(:,3) + planePosition(3);

% Draw the plane
tmpColor = ones(1, size(x, 2));
p = feval(@patch, x, y, z, tmpColor, 'edgeColor', 'none', patchArgs{:});
view(3);
