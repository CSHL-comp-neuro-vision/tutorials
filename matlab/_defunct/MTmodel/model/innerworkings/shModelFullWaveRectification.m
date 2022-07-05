% [pop, ind, res] = shmodv1simplerect(pop, ind, pars, resdirs)

function varargout = shmodv1simplerect(varargin)

% PARSE INPUTS
pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
if nargin > 3
    resdirs = varargin{4};
end

% OUTPUT
pop = pop.^2;
pop = pop * pars.scaleFactors.v1FullWaveRectified;
ind = ind;
varargout{1} = pop;
varargout{2} = ind;
if nargin > 3
    res = pop * pinv(shQwts(pars.v1PopulationDirections))' * shQwts(resdirs)';
    varargout{3} = res;
end