
function varargout = shModelV1Normalization(varargin)

% extract the necessary parameters
pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
 popdirs = pars.v1PopulationDirections;          % directions (in 3d Fourier space) of the v1 filters
 v1normtype = pars.v1NormalizationType;
if nargin > 3
    resdirs = varargin{4};
end

% call a the normalization routine specified by pars.v1normtype.
if strcmp(v1normtype, 'off')

    pop = pop;
    ind = ind;
    nume = pop;
    deno = 1;
    
    varargout{1} = pop;
    varargout{2} = ind;
    varargout{3} = nume;
    varargout{4} = deno;
    
    if nargin > 3
        res = pop * pinv(shQwts(popdirs))' * shQwts(resdirs)';     
        varargout{5} = res;
    end

elseif strcmp(v1normtype, 'global')        % OLD, UNIFORM NORMALIZATION
    
    if nargin < 4
        [pop, ind, nume, deno] = shModelV1Normalization_Global(pop, ind, pars);
        varargout{1} = pop;
        varargout{2} = ind;
        varargout{3} = nume;
        varargout{4} = deno;
    else
        [pop, ind, nume, deno, res, resnume, resdeno] = shModelV1Normalization_Global(pop, ind, pars, resdirs);
        varargout{1} = pop;
        varargout{2} = ind;
        varargout{3} = nume;
        varargout{4} = deno;
        varargout{5} = res;
        varargout{6} = resnume;
        varargout{7} = resdeno;
    end
        
elseif strcmp(v1normtype, 'tuned')                       % NEW, TUNED NORMALIZATION
    if nargin < 4
        [pop, ind, nume, deno] = shModelV1Normalization_Tuned(pop, ind, pars);
        varargout{1} = pop;
        varargout{2} = ind;
        varargout{3} = nume;
        varargout{4} = deno;
    else
        [pop, ind, nume, deno, res, resnume, resdeno] = shModelV1Normalization_Tuned(pop, ind, pars, resdirs);
        varargout{1} = pop;
        varargout{2} = ind;
        varargout{3} = nume;
        varargout{4} = deno;
        varargout{5} = res;
        varargout{6} = resnume;
        varargout{7} = resdeno;
    end
    
    % elseif strcmp(v1normtype, 'tuned2')
    %     [res, ind, pop, resnorm, Nnrmin, Nnrmout] = shModelV1Normalization_Tuned2(Qpop, Qind, Qres, pars, resdirs);
    %     popnorm = Nnrmin;
end
