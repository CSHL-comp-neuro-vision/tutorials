% [res, ind, pop, resnorm, popnorm, Qres, Qpop] = shModelV1Normalization_Tuned1(Qpop, Qind, Qres, pars, extradirs)
% do "type 1" tuned normalization on V1 neurons
 
function varargout = shModelV1Normalization_Tuned1(varargin)

pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
 v1sigma = pars.v1C50;         % semi-saturation constant of normalization signal
 popdirs = pars.v1PopulationDirections;           % directions (in 3d Fourier space) of the v1 filters
 scales = pars.nScales;           % number of spatial scales
 normstrength = pars.scaleFactors.v1NormalizationStrength;
if nargin > 3
    resdirs = varargin{4};
end

% get the filters ready. We have to reshape them into 3D filters.
xfilt = pars.v1NormalizationSpatialFilter;
tfilt = pars.v1NormalizationTemporalFilter;
trimmer = [length(xfilt), length(xfilt), length(tfilt)];
trimmer = trimmer - 1;
trimmer = trimmer./2;
deno = shGaussianBlur(pop, ind, xfilt, tfilt);

% now get normalizing
[pop, ind] = shTrim(pop, ind, trimmer);
nume = pop * pars.scaleFactors.v1Complex;
normwts = shModelV1Normalization_TunedWts(popdirs, popdirs)';
normwts = normwts*pars.scaleFactors.v1NormalizationPopulationK;
pop = nume./(normstrength.*deno*normwts + v1sigma.^2);

varargout{1} = pop;
varargout{2} = ind;
varargout{3} = nume;
varargout{4} = deno;

if nargin > 3
    normwts = shModelV1Normalization_TunedWts(resdirs, popdirs)';
    normwts = normwts*pars.scaleFactors.v1NormalizationPopulationK;
    
    resnume = nume * pinv(shQwts(popdirs))' * shQwts(resdirs)';
    resdeno = deno*normwts;

    res = resnume./(normstrength.*resdeno + v1sigma.^2);
    
    varargout{5} = res;
    varargout{6} = resnume;
    varargout{7} = resdeno;
end
