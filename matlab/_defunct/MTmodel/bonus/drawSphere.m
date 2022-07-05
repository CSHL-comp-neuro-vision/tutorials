% p = drawSphere(spherePosition, sphereRadius, nAngles, patchArgs)
%
% Draw a sphere on the current axes.
%
% Required arguments:
% sphereRadius          The radius of the sphere
%
% Optional arguments:
% spherePosition        The center point of the sphere. DEFAULT = the
%                       center point of the current axes.
% nAngles               The number of angles, both azimuthal and
%                       elevational, into which the sphere is divided.
%                       DEFAULT = 20
% patchArgs             Any number of arguments that are passed along to
%                       the patch command that actually draws the sphere, 
%                       e.g. 'facecolor', etc. DEFAULT = 'edgecolor',
%                       'none'

function p = drawSphere(varargin)

nAngles = 'default';
patchArgs = 'default';

sphereRadius = varargin{1};
if nargin >= 2;         spherePosition = varargin{2};   end
if nargin >= 3;         nAngles = varargin{3};          end
if nargin >= 4;         patchArgs = varargin(4:end);    end

if strcmp(spherePosition, 'default');
    xLim = get(gca, 'xlim');
    yLim = get(gca, 'ylim');
    zLim = get(gca, 'zlim');
    clear spherePosition
    spherePosition(1) = yLim(1) + (yLim(2)-yLim(1))/2;
    spherePosition(2) = xLim(1) + (xLim(2)-xLim(1))/2;
    spherePosition(3) = zLim(1) + (zLim(2)-zLim(1))/2;
end 
if strcmp(nAngles, 'default');          nAngles = 20;       end
if strcmp(patchArgs, 'default');        patchArgs = {};     end

phis = linspace(-pi/2, pi/2, nAngles+1);
thetas = linspace(0, 2.*pi, nAngles+1);
y = zeros(4, nAngles.^2);
x = y;
z = y;
rvec = [sphereRadius; sphereRadius; sphereRadius; sphereRadius];
for iTheta = 1:nAngles
    for jPhi = 1:nAngles
        az = [thetas(iTheta); thetas(iTheta+1); thetas(iTheta+1); thetas(iTheta)];
        el = [phis(jPhi); phis(jPhi); phis(jPhi+1); phis(jPhi+1)];
        yxz = sphere2rec([az, el, rvec]);
  
        
       
        pt = (iTheta-1).*nAngles+jPhi;
        y(:,pt) = yxz(:,1) + spherePosition(1);
        x(:,pt) = yxz(:,2) + spherePosition(2);
        z(:,pt) = yxz(:,3) + spherePosition(3);
    end
end

c = ones(nAngles.^2, 1);
p = feval(@patch, x, y, z, c', 'edgecolor', 'none', patchArgs{:});
