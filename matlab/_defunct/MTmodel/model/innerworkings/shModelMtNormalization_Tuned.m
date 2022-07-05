% [res, ind, pop, resnorm, popnorm, Qres, Qpop] = shModelV1Normalization_Tuned1(Qpop, Qind, Qres, pars, extradirs)
% do "type 1" tuned normalization on V1 neurons
 
function varargout = shModelV1Normalization_Tuned1(varargin)

pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
 mtsigma = pars.mtC50;                    % semi-saturation constant of normalization signal
 popvels = pars.mtPopulationVelocities;                     % directions (in 3d Fourier space) of the v1 filters
 scales = pars.nScales;                      % number of spatial scales
 normstrength = pars.scaleFactors.mtNormalizationStrength;        % strength of normalization signal
if nargin > 3
    res = varargin{4};
    resvels = varargin{5};
end
 
% get the filters ready. We have to reshape them into 3D filters.
xfilt = pars.mtNormalizationSpatialFilter;
tfilt = pars.mtNormalizationSpatialFilter;
trimmer = [length(xfilt), length(xfilt), length(tfilt)];
trimmer = trimmer - 1;
trimmer = trimmer./2;
deno = shGaussianBlur(pop, ind, xfilt, tfilt);


% now get normalizing
if nargin > 3
    resnume = shTrim(res, ind, trimmer);
    normwts = ones(size(res, 2), size(deno, 2))';
    resdeno = deno*normwts;
    
    res = pars.scaleFactors.mtPattern.*resnume./(normstrength.*resdeno + mtsigma.^2);
    
    
    varargout{5} = res;
    varargout{6} = resnume;
    varargout{7} = resdeno;
end

[pop, ind] = shTrim(pop, ind, trimmer);
nume = pop;
normwts = ones(size(pop, 2), size(deno, 2))';
pop = pars.scaleFactors.mtPattern.*pop./(normstrength.*deno*normwts' + mtsigma.^2);

varargout{1} = pop;
varargout{2} = ind;
varargout{3} = nume;
varargout{4} = deno;
