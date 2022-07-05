% [pop, ind, res] = shmodv1simplerect(pop, ind, pars, resdirs)

function varargout = shModelHalfWaveRectification(varargin)

% PARSE INPUTS
pop = varargin{1};
ind = varargin{2};
pars = varargin{3};

if nargin > 3
    res = varargin{4};
    res = res + pars.mtAlpha;
    res = (res>0).*res.^pars.mtExponent;
end

pop = pop + pars.mtAlpha;
pop = (pop>0).*pop.^pars.mtExponent;

ind = ind;

varargout{1} = pop.*pars.scaleFactors.mtHalfWaveRectification;
varargout{2} = ind;
if nargin > 3
    varargout{3} = res.*pars.scaleFactors.mtHalfWaveRectification;
end
