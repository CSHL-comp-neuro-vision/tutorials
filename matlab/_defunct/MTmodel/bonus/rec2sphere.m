% sphericalPoints = rec2sphere(rectangularPoints)     
%
% Transform [y, x, t] coordinates to [az, el, radius]
%
% Required arguments:
% rectangularPoints     the points in rectangular coordinates you want 
%                       transformed into spherical coordinates. Each row
%                       contains a different point in [Y X T] coordinates.
%
% Output:
% sphericalPoints       the transformed points. Each row of sphericalPoints
%                       corresponds to the same row of rectangularPoints.
%                       The first column specifies the azimuthal angle in
%                       radians with 0 = right. The second column specifies
%                       the elevation angle in radians from -pi/2 to pi/2,
%                       with 0 lying on the XY plane. The third column
%                       specifies the radius.

function sphericalPoints = rec2sphere(rectangularPoints)

sphericalPoints = zeros(size(rectangularPoints));

sphericalPoints(:,3) = sqrt(sum(rectangularPoints.^2,2));
sphericalPoints(:,1) = atan3(rectangularPoints(:,1), rectangularPoints(:,2));
sphericalPoints(:,2) = atan2(rectangularPoints(:,3), sqrt(sum(rectangularPoints(:,1:2).^2, 2)));