% crossAxes(varargin)   Put good old fashioned axes through the center of your plot.
%
% Sick of those awful box axes MATLAB is so in love with? Crossaxes is your
% solution. Feed it any PLOT arguments, such as 'linewidth' etc., to get
% them looking the way you want them.
%
% The varargin contains any valid arguments for PLOT commands, e.g.
% 'linestyle', 'color', etc.
%
% The axes are drawn at the center of the current axes, so if you want to
% change where they are placed, change the limits of the axes before
% calling crossAxes.


function crossAxes(varargin)

xl = get(gca, 'xlim');
yl = get(gca, 'ylim');
zl = get(gca, 'zlim');

xp = sum(xl)./2;
yp = sum(yl)./2;
zp = sum(zl)./2;

hold on

feval(@plot3, ...
    linspace(xl(1),xl(2)), ...
    yp*ones(1,100), ...
    zp*ones(1,100), ...
    varargin{:});
        
feval(@plot3, ...
    xp*ones(1,100), ...
    linspace(yl(1),yl(2)), ...
    zp*ones(1,100), ...
    varargin{:});

feval(@plot3, ...
    xp*ones(1,100), ...
    yp*ones(1,100), ...
    linspace(zl(1),zl(2)), ...
    varargin{:});

hold off
