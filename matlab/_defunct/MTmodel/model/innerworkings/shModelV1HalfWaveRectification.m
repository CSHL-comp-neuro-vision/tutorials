% [pop, ind, res] = shmodv1simplerect(pop, ind, pars, resdirs)

function varargout = shmodv1simplerect(varargin)

% PARSE INPUTS
pop = varargin{1};
ind = varargin{2};
pars = varargin{3};

if nargin > 3
    resdirs = varargin{4};
    res = shSwts(resdirs) * pinv(shSwts(pars.v1dirs)) * pop;
    res = (res>0).*res.^2;
end
pop = (pop>0).*pop.^2;

ind = ind;
varargout{1} = pop;
varargout{2} = ind;
if nargin > 3 
    varargout{3} = res;
end