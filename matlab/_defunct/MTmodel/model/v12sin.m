% function preferredSin = v12sin(v1Neurons)  
%
% Get the paramters of the drifting grating preferred by given V1 neurons.
%
% Required argument:
% v1Neurons             a matrix containing the parameters of the neurons
%                       whose preferred gratings you wish to find. Each row
%                       of v1Neurons specifies a neuron in the standard
%                       way: the first column gives the preferred
%                       direction, and the second column gives the
%                       preferred ratio of temporal frequency to spatial
%                       frequency in cycles/frame and cycles/pixel,
%                       respectively.
%
% Output:
% preferredSin          a matrix containing the parameters of the gratings
%                       preferred by the neurons given in v1Neurons. Each
%                       row gives the parameters for the grating preferred
%                       by the neuron given by the corresponding row of
%                       v1Neurons. The first column gives the preferred
%                       direction in radians with 0 = right. The second
%                       column gives the preferred spatial frequency in
%                       cycles/pixel. The third column gives the preferred
%                       temporal frequency in cycles/frame.

function preferredSin = v12sin(v1Neurons)

k = .2173.*ones(size(v1Neurons,1),1);
preferredSin = zeros(size(v1Neurons,1), 3);

% now make the conversion

preferredSin(:,1) = v1Neurons(:,1);
preferredSin(:,2) = k.*cos(atan3(v1Neurons(:,2), ones(size(v1Neurons,1), 1)));
preferredSin(:,3) = k.*sin(atan3(v1Neurons(:,2), ones(size(v1Neurons,1), 1)));