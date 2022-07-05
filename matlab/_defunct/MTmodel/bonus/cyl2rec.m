% recPoints = sphere2rec(cylPoints)     Transform [az, h, r] coordinates to [y, x, t]
%
% cylPoints is a matrix, each row of which describes a point in 3-space in
% [azimuthal angle, height, radius] coordinates. The azimuthal angle is in
% radians, with 0 = right.
% 
% recPoints is a matrix each row of which describes the same point as the
% corresponding row of cylPoints, but in [Y, X, T] coordinates.

function recPoints = sphere2rec(cylPoints)

if size(cylPoints,2) == 2;
    cylPoints = [cylPoints, ones(size(cylPoints,1), 1)];
end

recPoints = zeros(size(cylPoints));

recPoints(:,1) = cylPoints(:,3).*sin(cylPoints(:,1));
recPoints(:,2) = cylPoints(:,3).*cos(cylPoints(:,1));
recPoints(:,3) = cylPoints(:,2);