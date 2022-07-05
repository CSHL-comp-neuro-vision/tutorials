% labelCrossAxes(xLabel, yLabel, zLabel, textArgs)    Label crossAxes
%
% Required arguments:
% xLabel            the label for the x-axis
% yLabel            the label for the y-axis
% zLabel            the label for the z-axis
%
% Optional arguments:
% labelPosition     where the labels should be placed, in units normalized
%                   to the extent of the current axes. DEFAULT = 1.1
% textArgs          any number of arguments that are passed on to the TEXT
%                   command that writes the labels on the axes.

function labelCrossAxes(varargin)

labelPosition = 'default';
textArgs = 'default';

xLabel = varargin{1};
yLabel = varargin{2};
zLabel = varargin{3};
if nargin >= 4;         labelPosition = varargin{4};            end
if nargin >= 5;         textArgs = varargin(5:end);             end

if strcmp(labelPosition, 'default');        labelPosition = [1.1 1.1 1.1];      end
if strcmp(textArgs, 'default');             textArgs = {};                      end


% DONE PARSING ARGUMENTS
xLim = get(gca, 'xlim');
yLim = get(gca, 'ylim');
zLim = get(gca, 'zlim');

xPosition = sum(xLim)./2;
yPosition = sum(yLim)./2;
zPosition = sum(zLim)./2;

feval(@text, labelPosition(1)*xLim(2), yPosition, zPosition, xLabel, textArgs{:})
feval(@text, xPosition, labelPosition(2)*yLim(2), zPosition, yLabel, textArgs{:})
feval(@text, xPosition, yPosition, labelPosition(3)*zLim(2), zLabel, textArgs{:})
