% res = shGetNeuron(shMatrix, ind, neuronsToGet, scaleToGet, neuronPosition) 
%
% Extract the response of neuron(s) at one spatial position from shModel outputs.
%
% Required arguments:
% shMatrix          output from the shModel function. Either the POP or RES
%                   outputs are acceptable.
% ind               the IND output from the shModel function.
% 
% Optional arguments:
% neuronsToGet      a vector whose elements point to which neurons from
%                   shMatrix you want to get. For instance, if shMatrix
%                   contains the responses of MT neurons tuned for [0 1],
%                   [0 4], [pi 2], neuronsToGet could be [2, 3]. In that
%                   case the output would contain two rows, the first row
%                   containing the responses of a neuron tuned to [0 4],
%                   the second row containing the responses of a neuron
%                   tuned to [pi 2]. DEFAULT = you the responses of one
%                   neuron from each set of similarly tuned neurons.
% scaleToGet        not currently implemented. DEFAULT = 1.
% neuronPosition    the spatial position of the neurons whose responses you
%                   want, in [Y X] coordinates. [0 0] is a point closest to
%                   the center of the stimulus. DEFAULT = [0 0].

function res = shGetNeuron(varargin)

neuronsToGet = 'default';
scaleToGet = 'default';
neuronPosition = 'default';

shMatrix = varargin{1};
ind = varargin{2};
if nargin >= 3;         neuronsToGet = varargin{3};         end
if nargin >= 4;         scaleToGet = varargin{4};           end
if nargin >= 5;         neuronPosition = varargin{5};       end

if strcmp(neuronsToGet, 'default');     neuronsToGet = [1:size(shMatrix, 2)];       end;
if strcmp(scaleToGet, 'default');       scaleToGet = 1;                             end;
if strcmp(neuronPosition, 'default');   neuronPosition = [0 0];                     end;

% DONE PARSING INPUTS

for i = 1:length(neuronsToGet)
    tmp = shGetSubPop(shMatrix, ind, neuronsToGet(i), scaleToGet);
    cpy = floor(size(tmp,1)./2) + 1;
    cpx = floor(size(tmp,2)./2) + 1;

    res(i,:) = squeeze(tmp(cpy + neuronPosition(1), cpx + neuronPosition(2), :))';
end
