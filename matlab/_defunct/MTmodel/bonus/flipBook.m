% flipBook(matrixToFlip, displayRange, pauseLength)       
%
% Show the values of a 3D matrix as a movie.
%
% Required arguments: 
% matrixToFlip          the 3D matrix to be shown as a movie
%
% Optional arguments:
% displayRange          the min and max values to show. DEFAULT = the min
%                       and max values of matrixToFlip.
% pauseLength           the time to pause before displaying the next frame
%                       in seconds. DEFAULT = 0.

function flipBook(varargin)

displayRange = 'default';
pauseLength = 'default';

matrixToFlip = varargin{1};
if nargin >= 2;         displayRange = varargin{2};         end
if nargin >= 3;         pauseLength = varargin{3};          end

if strcmp(displayRange, 'default')
    displayRange = [min2(matrixToFlip), max2(matrixToFlip)];
end
if strcmp(pauseLength, 'default');      pauseLength = 0;        end

% DONE PARSING ARGUMENTS. NOW FLIP THAT BOOK!
for i = 1:size(matrixToFlip, 3);
    showIm(matrixToFlip(:,:,i), displayRange, 'auto', 0);
    drawnow
    pause(pauseLength);
end
