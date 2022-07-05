% res = shGetSubPop(shMatrix, ind, neuronsToGet, scaleToGet)
%
% Extract the responses of all the neurons with similar tuning at different
% spatial positions from the output of shModel.
%
% Extract the response of neuron(s) at one spatial position from shModel outputs.
%
% Required arguments:
% shMatrix          output from the shModel function. Either the POP or RES
%                   outputs are acceptable.
% ind               the IND output from the shModel function.
% 
% Optional arguments:
% neuronsToGet      a number that points to which neuron from
%                   shMatrix you want to get. For instance, if shMatrix
%                   contains the responses of MT neurons tuned for [0 1],
%                   [0 4], [pi 2], neuronsToGet could be 2. In that
%                   case the output would be a 3D matrix whose entries
%                   contained the responses of all the neurons tuned to [0
%                   4], but with RF centered on different spatial
%                   positions. DEFAULT = 1.
% scaleToGet        not currently implemented. DEFAULT = 1.

function res = shGetSubPop(varargin)

neuronsToGet = 'default';
scaleToGet = 'default';

shMatrix = varargin{1};
ind = varargin{2};
if nargin >= 3;         neuronsToGet = varargin{3};             end;
if nargin >= 4;         scaleToGet = varargin{4};               end;

if strcmp(neuronsToGet, 'default');         neuronsToGet = 1;       end;
if strcmp(scaleToGet, 'default');           scaleToGet = 1;         end;    

scaleback = 1;

shMatrix = shMatrix(ind(scaleToGet,1)+1:ind(scaleToGet+1,1), neuronsToGet);
res = reshape(shMatrix, [ind(scaleToGet+1, 2), ind(scaleToGet+1, 3), ind(scaleToGet+1, 4), 1]);

if scaleback == 1;
    while scaleToGet > 1
        x = size(res, 1);
        x1 = reshape(linspace(0, 1, x), [x 1 1]);
        x2 = reshape(linspace(0, 1, 2*x), [2*x 1 1]);
        y = size(res, 2);
        y1 = reshape(linspace(0, 1, y), [1 y 1]);
        y2 = reshape(linspace(0, 1, 2*y), [1 2*y 1]);
        t = size(res, 3);
        t1 = reshape(linspace(0, 1, t), [1 1 t]);
        t2 = reshape(linspace(0, 1, 2*t), [1 1 2*t]);

        [X1, Y1, T1] = ndgrid(x1, y1, t1);
        [X2, Y2, T2] = ndgrid(x2, y2, t2);

        res = interpn(X1, Y1, T1, res, X2, Y2, T2);

        scaleToGet = scaleToGet-1;
    end
end
