
function varargout = shmodmtpostfilt(varargin)

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
if pars.mtSpatialPoolingBeforeThreshold == 0
    f = pars.mtSpatialPoolingFilter;
    
    if nargin > 3
        res = shGaussianBlur(res, ind, f);
        varargout{3} = res;
    end
    
    [pop, ind] = shGaussianBlur(pop, ind, f);
    varargout{1} = pop;
    varargout{2} = ind;
end


varargout{1} = pop;
varargout{2} = ind;
if nargin > 3
    varargout{3} = res;
end
