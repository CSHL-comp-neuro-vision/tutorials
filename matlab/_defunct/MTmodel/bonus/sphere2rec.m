% rectangularPoints = sphere2rec(sphericalPoints)     
%
% Transform [az, el, radius] coordinates to [y, x, t]
%
% Required arguments:
% sphericalPoints       The points in spherical coordinates you want
%                       transformed into rectangular coordinates. Each row
%                       contains a different point. The first column
%                       specifies the azimuthal angle in radians with 
%                       0 = right. The second column specifies
%                       the elevation angle in radians from -pi/2 to pi/2,
%                       with 0 lying on the XY plane. The third column
%                       specifies the radius.
%
% Output:
% rectangularPoints     the transformed points Each row contains a different
%                       point in [Y X T] coordinates.


function rectangularPoints = sphere2rec(sphericalPoints)

if size(sphericalPoints,2) == 2;
    sphericalPoints = [sphericalPoints, ones(size(sphericalPoints,1), 1)];
end

rectangularPoints = zeros(size(sphericalPoints));

rectangularPoints(:,1) = sphericalPoints(:,3).*cos(sphericalPoints(:,2)).*sin(sphericalPoints(:,1));
rectangularPoints(:,2) = sphericalPoints(:,3).*cos(sphericalPoints(:,2)).*cos(sphericalPoints(:,1));
rectangularPoints(:,3) = sphericalPoints(:,3).*sin(sphericalPoints(:,2));