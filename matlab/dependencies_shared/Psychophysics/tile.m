function tile(m,n,monitor,menuBarFlag)
%tile([m],[n],[monitor],[menuBarFlag])
%
%Tiles matlab figures in the monitor in rows starting at the top-right.
%
%Inputs:
%   m  #of rows of figures (default is 3)
%   n  #of columns (default is 3)
%   monitor monitor number (1 is desktop, 2 is extended, default is
%       highest available)
%   menuBarFlag determines if the menu bar is shown.  Default is 1 (on)

%Written 3/30/12 gmb (after getting frustrated by Matlab's positioning of
%figures outside the monitor by default)



monitorPos = get(0,'MonitorPositions');

nMonitors = size(monitorPos,1);

if ~exist('monitor','var')
    monitor= nMonitors;
end

if ~exist('menuBarFlag')
    menuBarFlag = 1;
end

if monitor>nMonitors
    error(sprintf('No monitor #%d available',monitor));
end

sz = monitorPos(monitor,[3,4])-monitorPos(monitor,[1,2])+1;

if ~exist('m','var')
    m = 3;
end

if ~exist('n','var')
    n = 3;
end

x = 0;

if monitor==1
    x0=-2;
    y0 = sz(2)+2;
    dx = round((sz(1)+2)/n);
dy = round((sz(2)-38)/m);
else
    x0 = monitorPos(2,1)-2;
    y0 = monitorPos(1,4)-monitorPos(2,2)+2;
    dx = round((sz(1)+2)/n);
dy = round(sz(2)/m);
end
x =x0;
y =y0;


figs = sort(get(0,'Children'))';
set(0,'Units','pixels');

for i=figs
    set(i,'units','pixels');
    set(i,'OuterPosition',[x,y-dy,dx,dy]);
    if menuBarFlag
        set(i,'MenuBar','figure');
    else
        set(i,'MenuBar','none');
    end
    x = x+dx;
    if x>=x0+sz(1)-5
        x = x0;
        y = y-dy;
    end
end

