% [res, pop] = shmodv1comp(N, pars);      get model complex cell responses
% N is a matrix containing the firing rates of model v1 neurons, following
% squaring and normalization. pop is a matrix containing the four-dimensional
% responses of model complex cells.

function varargout = shModelV1Blur(varargin)

pop = varargin{1};
ind = varargin{2};
pars = varargin{3};
 filt = pars.v1ComplexFilter;
 popdirs = pars.v1PopulationDirections;
 nScales = pars.nScales;
if nargin > 3
    resdirs = varargin{4};
end

% 4. convolve each neuron's response with the blurring filter to get local
% motion energy: i.e., model complex cell responses    
[pop, ind] = shGaussianBlur(pop, ind, filt); 
pop = pop*pars.scaleFactors.v1Blur;

% 5. Steer population responses to get responses for neurons specified in
% resdirs.
varargout{1} = pop;
varargout{2} = ind;
if nargin > 3
    res = pop * pinv(shQwts(popdirs))' * shQwts(resdirs)';
    varargout{3} = res;
end
