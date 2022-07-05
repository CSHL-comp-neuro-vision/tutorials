% thisScale = shGetScale(shMatrix, ind, scaleToGet)   
%
% Extract all the neuronal responses at a particular scale from shMatrix,
% which is the output of a call to shModel. The responses will be in 4D
% format with the coordinate system [Y X T N].
%
% Required arguments:
% shMatrix          the matrix from which you want to extract values. This
%                   can be either the POP or RES outputs of shModel.
% ind               the IND output of the same call to shModel
% 
% Optional arguments:
% scaleToGet        the scale for which you want the values. DEFAULT = 1
%
% Output
% thisScale         the extracted values in [Y X T N] format. Y and X refer
%                   to the spatial center of the neuron's receptive field.
%                   T is the time. N is the tuning of the neuron.


function thisScale = shGetScale(shMatrix, ind, scaleToGet)

if exist('scaleToGet') ~= 1
    scaleToGet = 1;
end

shMatrix = shMatrix(:, ind(scaleToGet)+1:ind(scaleToGet+1))';
thisScale = reshape(shMatrix, [ind(scaleToGet+1, 2), ind(scaleToGet+1, 3), ind(scaleToGet+1, 4), size(shMatrix, 2)]);
