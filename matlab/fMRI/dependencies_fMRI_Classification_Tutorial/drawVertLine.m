function drawVertLine(locs,color,style)

%-------------------------------------------
%
% drawVertLine(locs,color,style)
%
% simple function for drawing vertical lines
% in the current axis
%
% required input:
% locs -- a scalar or vector, where to draw the lines
%
% optional inputs:
% color -- line color (default = 'k')
% style -- line style (default = '-')
%
% freeman, 10-25-2009
%-------------------------------------------

if ~exist('color','var');
    color = 'k';
end

if ~exist('style','var');
    style = '-';
end

% get the current axis so we know the ylim
tmpProp = get(gca);
ylim = tmpProp.YLim;

hold on
for iLine = 1:length(locs)
	h = line([locs(iLine) locs(iLine)],[ylim(1) ylim(2)]);
    set(h,'Color',color)
    set(h,'LineStyle',style)
end
hold off