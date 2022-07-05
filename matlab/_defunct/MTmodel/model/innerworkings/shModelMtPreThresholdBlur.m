% [res, pop] = shModelMtLinear(C, pars, resvels)   get linear MT responses from complex cell responses
% C is returned by shmodv1comp. pop is the response of MT neurons before
% squaring and normalization. pop contains the responses of neurons
% preferring velocities specified by pars.popvels. The user can get the
% responses of additional neurons by supplying resvels, which must be
% an mx2 matrix specifying the preferred velocities of m extra neurons. If
% resvels is supplied, res contains the firing rates of the neurons
% specified by resvels; if it is not supplied, res is a copy of pop.

function varargout = shModelMtLinear(varargin)

pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
v1normtype = pars.v1NormalizationType;
v1dirs = pars.v1PopulationDirections;
popvels = pars.mtPopulationVelocities;
nscales = pars.nScales;
if nargin > 3
    res = varargin{4};
end


% if PARS calls for spatial pooling before the pattern computation, do it now.
if pars.mtSpatialPoolingBeforeThreshold == 1
    f = pars.mtSpatialPoolingFilter;

    if nargin > 3
        res = shGaussianBlur(res, ind, f);
        varargout{3} = res;
    end

    [pop, ind] = shGaussianBlur(pop, ind, f);
end

varargout{1} = pop;
varargout{2} = ind;
if nargin > 3
    varargout{3} = res;
end
