% [xDirection, yResponse] = shTunePlaidOriOri(pars, neuron, stageName, reso)  
%
% Generate an ori/ori surface for a model neuron. A neuron is shown a
% series of plaids whose components' directions of motion both vary
% independently. The reponse to each plaid is stored in a matrix and viewed
% graphically as a surface with the directions of the two components being
% the independent variables.
%
% This code also computes the respone of the neuron to individual gratings
% at half contrast. These are stored in the last row and the first column.
% 
% Output:
% xDirection        a vector of the directions that each component can take
% yResponse         a matrix containing the responses of the neuron to
%                   each stimulus condition
%
% Required arguments:
% pars              a parameters structure
% neuron            parameters of the neuron you want to model
% stageName         name of the stage you want to model. Choices:
%                   'v1Complex' and 'mtPattern'
%
% Optional arguments:
% nDataPoints       the number of different stimulus conditions. 
%                   DEFAULT = 7
% gratingSf         the spatial frequency of the plaid's grating components,
%                   in cycles/pixel. DEFAULT = the neuron's preferred 
%                   spatial frequency.
% gratingTf         the temporal frequency of the plaid's grating components,
%                   in cycles/frame. DEFAULT = the neuron's preferred temporal
%                   frequency.
% plaidAngle        the angle between the plaid's component gratings, in
%                   radians. DEFAULT = (2/3)*pi  (120 degrees).
% plaidContrast     the contrast of the plaid. DEFAULT = pars.v1C50.


function [xDirection, yResponse] = shTuneGratingDirectionori(varargin)

% The following arguments are optional and are by default 'default'.
nDataPoints = 'default';
gratingSf = 'default';
gratingTf = 'default';
plaidAngle = 'default';
plaidContrast = 'default';

% parse the varargin
                    pars = varargin{1};
                    neuron = varargin{2};
                    stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};              end
if nargin >= 5;     gratingSf = varargin{5};                end
if nargin >= 6;     gratingTf = varargin{6};                end
if nargin >= 7;     plaidAngle = varargin{7};               end
if nargin >= 8;     plaidContrast = varargin{8};            end

if strcmp(stageName(1:2), 'v1');
    preferredGrating = v12sin(neuron);
else
    preferredGrating = mt2sin(neuron);
end

% Assign default values
if strcmp(nDataPoints, 'default');          nDataPoints = 11;                       end
if strcmp(gratingSf, 'default');            gratingSf = preferredGrating(2);        end
if strcmp(gratingTf, 'default');            gratingTf = preferredGrating(3);        end
if strcmp(plaidAngle, 'default');           plaidAngle = (2/3)*pi;                  end
if strcmp(plaidContrast, 'default');        plaidContrast = pars.v1C50;             end

% done parsing arguments. Now get on with it.
stimSz = shGetDims(pars, stageName, [1 1 31]);
xDirection = linspace(0, 2*pi, nDataPoints+1);
xDirection = xDirection(1:end-1);
xDirection = [-1, xDirection];
yResponse = zeros(length(xDirection), length(xDirection));
for i = 1:length(xDirection)
    for j = i:length(xDirection)
        
        % Generate the stimulus
        if i == 1;
            firstGrating = zeros(stimSz);
        else
            firstGrating = mkSin(stimSz, xDirection(i), gratingSf, ...
                                 gratingTf, plaidContrast/2);
        end
        if j == 1;
            secondGrating = zeros(stimSz);
        else
            secondGrating = mkSin(stimSz, xDirection(j), gratingSf, ...
                                 gratingTf, plaidContrast/2);
        end
        thisPlaid = firstGrating + secondGrating;
 
        % compute the neuron's response to this stimulus
        [pop, ind, res] = shModel(thisPlaid, pars, stageName, neuron);
        if strcmp(stageName, 'v1lin')
            res = sqrt(res.^2);
        end
        res = mean(shGetNeuron(res, ind));
        yResponse(end-j+1, i) = res;
        yResponse(end-i+1, j) = res;
        
        % show the results thus far
        showIm(yResponse);
        drawnow
    end
end