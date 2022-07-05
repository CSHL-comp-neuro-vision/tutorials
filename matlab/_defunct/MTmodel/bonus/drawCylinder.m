% p = drawCylinder(cylPosition, cylHeight, cylRadius, 
%                  azAngle, elAngle, nSides, varargin)
%
% Draw a cylinder on the current axes.
%
% Required arguments:
% cylPosition   the position of the cylinder in [y, x, t] coordinates
% cylHeight     the height of the cylinder
% cylRadius     the radius of the cylinder
% azAngle       the azimuthal angle of the cylinder's orientation, in
%               radians with 0 = right.
% elAngle       the elevation angle of the cylinder's orientation, in
%               radians from -pi/2 to pi/2, with 0 = lying on the X-Y plane
% nSides        the number of sides to draw. The higher the number, the
%               rounder the cylinder's appearance.
% varargin      these can contain any argument pairs that you can pass to
%               the patch function, such as 'color', 'facealpha', etc.
%
% Output:
% 
%
% SEE ALSO: PATCH

function p = drawSphere(cylPosition, cylHeight, cylRadius, azAngle, elAngle, nSides, varargin)

azrot = [cos(azAngle) sin(azAngle) 0; -sin(azAngle) cos(azAngle) 0; 0 0 1];
elrot = [1 0 0; 0 -sin(elAngle) cos(elAngle); 0 cos(elAngle) sin(elAngle)]; 

theta = linspace(0, 2.*pi, nSides+1);
y = zeros(4, nSides);
x = y;
z = y;
rvec = [cylRadius; cylRadius; cylRadius; cylRadius];
for i = 1:nSides
        az = [theta(i); theta(i+1); theta(i+1); theta(i)];
        hvec = cylHeight./2.*[-1; -1; 1; 1];
        yxz = cyl2rec([az, hvec, rvec]);
        
        yxz = (azrot*elrot*yxz')';
        
        y(:,i) = yxz(:,1) + cylPosition(1);
        x(:,i) = yxz(:,2) + cylPosition(2);
        z(:,i) = yxz(:,3) + cylPosition(3);
end

tmpColor = ones(1, size(x, 2));
p = feval(@patch, x, y, z, tmpColor, 'edgecolor', 'none', varargin{:});
view(3);