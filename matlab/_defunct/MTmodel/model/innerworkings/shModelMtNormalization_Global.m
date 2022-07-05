% [res, ind, pop, resnorm, popnorm, Qres, C] = shModelV1Normalization_untuned(C, Cind, Qres, pars)
% do untuned normalization on V1 responses.

function varargout = shModelV1Normalization_Global(varargin)

% Extract parameters
pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
mtsigma = pars.mtsigma;                    % semi-saturation constant of normalization signal
popvels = pars.mtvels;                     % directions (in 3d Fourier space) of the v1 filters
scales = pars.scales;                      % number of spatial scales
normstrength = pars.mtnormstrength;        % strength of normalization signal
if nargin > 3
    res = varargin{4};
    resvels = varargin{5};
end


% get the normalization signal.
popnorm = zeros(scales, ind(2, end));
for s = 1:scales
    tmp = shGetScale(pop, ind, s);
    tmp = mean(mean(mean(tmp, 1), 2), 4);
    tmp = squeeze(tmp)';
    popnorm(s, 1:size(tmp, 2)) = tmp;
end

% normalize the population of responses
nume = zeros(size(pop));
deno = nume;
for s = 1:scales
    scale = shGetScale(pop, ind, s);
    scalenorm = popnorm(s, 1:ind(s+1, 4));    % get the normalization signal for the current s
    scalenorm = reshape(scalenorm, [1 1 size(scalenorm, 2) 1]);
    scalenorm = repmat(scalenorm, [size(scale, 1), size(scale, 2), 1, size(scale, 4)]);    % make it the right size

    scale = scale./(normstrength.*scalenorm + mtsigma);   % get this scale
    pop = shSetScale(pop, ind, scale, s);

    nume = shSetScale(nume, ind, scale, s);
    deno = shSetScale(deno, ind, scalenorm, s);

    if nargin > 3
        scale = shGetScale(res, ind, s);

        scalenorm = popnorm(s, 1:ind(s+1, 4));    % get the normalization signal for the current s
        scalenorm = reshape(scalenorm, [1 1 size(scalenorm, 2) 1]);
        scalenorm = repmat(scalenorm, [size(scale, 1), size(scale, 2), 1, size(scale, 4)]);    % make it the right size

        scale = scale./(normstrength.*scalenorm + mtsigma);
        res = shSetScale(res, ind, scale, s);
    end
end

% return stuff. Maybe do some interpolation to get extra responses.
varargout{1} = pop;
varargout{2} = ind;
varargout{3} = nume;
varargout{4} = deno;

if nargin > 3
    varargout{5} = res;
end
