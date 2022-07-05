function varargout = shModelMtLinear(varargin)

pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
popvels = pars.mtPopulationVelocities;
nscales = pars.nScales;
if nargin > 3
    resvels = varargin{4};
end


varargout{1} = pop*shMtWts(popvels, pars)'.*pars.scaleFactors.mtLinear;
varargout{2} = ind;
if nargin > 3
    varargout{3} = pop*shMtWts(resvels, pars)'.*pars.scaleFactors.mtLinear;
end
