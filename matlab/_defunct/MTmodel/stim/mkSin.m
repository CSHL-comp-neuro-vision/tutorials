% S =  mkSin(stimSz, sinDirection, sinSf, sinTf, sinContrast, sinPhase)    
% 
% mkSing returns a drifting grating. This grating will have a mean value of
% .5 and a contrast of sinContrast.
%
% Required arguments:
% stimSz            the size of the stimulus in [y, x, t] coordinates.
% sinDirection      the direction, in radians (0 = rightward), of motion
% sinSf             spatial frequency in units of cycles/pixel
% sinTf             temporal frequency in units of cycles/frame
%
% Optional arguments:
% sinContrast       grating contrast. DEFAULT = 1.
% sinPhase          initial phase in periods. DEFAULT = 0.

function [res] = mkSin(varargin)

% the following variables are optional and by default are 'default'
sinContrast = 'default';
sinPhase = 'default';

% parse varargin
                    stimSz = varargin{1};
                    sinDirection = varargin{2};
                    sinSf = abs(varargin{3});
                    sinTf = abs(varargin{4});
if nargin >= 5;     sinContrast = varargin{5};      end
if nargin >= 6;     sinPhase = varargin{6};         end

% assign default values
if strcmp(sinContrast, 'default');      sinContrast = 1;            end
if strcmp(sinPhase, 'default');         sinPhase = 0;               end

% Error message if arguments are incorrectly formatted
if length(stimSz) ~= 3
    error('stimSz must be a 3-vector');
end

% Make a coordinate system
y = [1:stimSz(1)] - (floor(stimSz(1)/2)+1);
x = [1:stimSz(2)] - (floor(stimSz(2)/2)+1);
t = [0:stimSz(3)-1];

[y, x, t] = ndgrid(x, y, t);
y = -y;

res = cos(2*pi*sinSf*cos(sinDirection)*x + ...
          2*pi*sinSf*sin(sinDirection)*y - ...
          2*pi*sinTf*t + ...
          2*pi*sinPhase);
          
res = sinContrast.*res./2 + .5;