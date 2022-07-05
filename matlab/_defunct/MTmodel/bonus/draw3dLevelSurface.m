% p = draw3dLevelSurface(f, levelToDraw, levelMode, color, alpha)
%
% Draw a 3D level surface of a 3D function.
%
% Required arguments:
% f             a matrix containing the values of a 3D function in [Y, X, T]
%               coordinates.
%
% Optional arguments:
% levelToDraw   The value of the function at which you want a level surface
%               drawn. DEFAULT = .5
% levelMode     if levelMode is 'norm', the levelToDraw level will be
%               relative to the maximum value of the function. Otherwise
%               the value will just be the absolute value of the function.
%               DEFAULT = 'norm'
% color         The color of the plotted surface.   DEFAULT = [0 0 1]
% alpha         The transparency of the surface, from 0 to 1. DEFAULT = 1


function p = draw3dLevelSurface(varargin);

levelToDraw = 'default';
levelMode = 'default';
color = 'default';
alpha = 'default';

                    f = varargin{1};
if nargin >= 2      levelToDraw = varargin{2};      end
if nargin >= 3      levelMode = varargin{3};        end
if nargin >= 4      color = varargin{4};            end
if nargin >= 5      alpha = varargin{5};            end

if strcmp(levelToDraw, 'default');      levelToDraw = .5;           end
if strcmp(levelMode, 'default');        levelMode = 'norm';         end
if strcmp(color, 'default');            color = [0 0 1];            end
if strcmp(alpha, 'default');            alpha = 1;                  end

%%%%%%%%%%%%%%%%%% END PARSING PARAMETERS

if strcmp(levelMode, 'norm')
    f = f./max2(f);
end

%%%%%%%%%%%%%%%%%% NOW GET PLOTTING

yl = get(gca, 'ylim');
xl = get(gca, 'xlim');
zl = get(gca, 'zlim');
y = linspace(yl(1), yl(2), size(f, 1));
x = linspace(xl(1), xl(2), size(f, 2));
z = linspace(zl(1), zl(2), size(f, 3));
[x,y,z] = meshgrid(x, y, z);

tmp = z;
for i = 1:size(y, 3);
    y(:,:,i) = flipud(y(:,:,i));
    z(:,:,i) = tmp(:,:,end-i+1);
end

p = patch(isosurface(x,y,z, f, levelToDraw));
isonormals(x,y,z,f,p)
set(p,'FaceColor', color,'EdgeColor','none', 'facealpha', alpha);
daspect([1 1 1]);
view(3);