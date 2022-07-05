function drawHorzLine(locs,color,style)

%-------------------------------------------
%
% drawHorzLine(locs,color,style)
%
% simple function for drawing a horizontal line
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

% get the current axis so we know the xlim
tmpProp = get(gca);
xlim = tmpProp.XLim;

hold on
for iLine = 1:length(locs)
	h = line([xlim(1) xlim(2)],[locs(iLine) locs(iLine)]);
    set(h,'Color',color)
    set(h,'LineStyle',style)
end
hold off
