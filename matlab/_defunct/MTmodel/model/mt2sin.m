% function preferredSin = mt2sin(mtNeurons)   
%
% Get the paramters for gratings preferred by given MT neurons.
%
% Required arguments:
% mtNeurons         an Nx2 matrix. Each row contains the parameters of a
%                   model MT neuron in standard form: the first column is
%                   the preffered direction in radians with 0 = right; the
%                   second column is the preferred speed in pixels/frame.
%
% Output:
% preferredSin      an Nx3 matrix. Each row contains the parameters of the
%                   sinusoid preferred by the MT neuron in the
%                   corresponding row of mtNeurons. The first column
%                   contains the preferred direction in radians with 0 =
%                   right; the second column contains the preferred spatial
%                   frequency in cycles/pixel; the third column contains
%                   the preferred temporal frequency in cycles/frame.

function preferredSin = mt2sin(mtNeurons)

r = .2173;
mtNeurons = [mtNeurons, ones(size(mtNeurons,1), 1)];

preferredSin(:,1) = mtNeurons(:,1);
preferredSin(:,2) = r.*cos(atan3(mtNeurons(:,2), ones(size(mtNeurons, 1), 1)));
preferredSin(:,3) = r.*sin(atan3(mtNeurons(:,2), ones(size(mtNeurons, 1), 1)));