% [pop, nrm] = shModelMtNormalization(linpop, pars)  get squared and normalized linpop responses
% linpop is the response returned by shModelMtLinear; pars is the SH parameters
% structure given by, e.g. shPars. pop contains the normalized
% responses; nrm contains the normalization signal.
% Extraneurons is a four-dimensional matrix containing the responses of
% neurons not in the the population

function varargout = shModelMtNormalization(varargin)

pop = varargin{1};
ind = varargin{2};
pars = varargin{3};

if nargin > 3
    res = varargin{4};
    resvels = varargin{5};
end

if strcmp(pars.mtNormalizationType, 'global')
    if nargin < 4
        [pop, ind, nume, deno] = shModelMtNormalization_Global(pop, ind, pars);
    else
        [pop, ind, nume, deno, res, resnume, resdeno] = shModelMtNormalization_Global(pop, ind, pars, res, resvels);
    end

elseif strcmp(pars.mtNormalizationType, 'tuned')
    if nargin < 4
        [pop, ind, nume, deno] = shModelMtNormalization_Tuned(pop, ind, pars);
    else
        [pop, ind, nume, deno, res, resnume, resdeno] = shModelMtNormalization_Tuned(pop, ind, pars, res, resvels);
    end

elseif strcmp(pars.mtNormalizationType, 'self')
    if nargin < 4
        [pop, ind, nume, deno] = shModelMtNormalization_Self(pop, ind, pars);
    else
        [pop, ind, nume, deno, res] = shModelMtNormalization_Self(pop, ind, pars, res, resvels);
    end
end


if nargin < 4
    varargout{1} = pop;
    varargout{2} = ind;
    varargout{3} = nume;
    varargout{4} = deno;
else
    varargout{1} = pop;
    varargout{2} = ind;
    varargout{3} = nume;
    varargout{4} = deno;
    varargout{5} = res;
    varargout{6} = resnume;
    varargout{7} = resdeno;
end
