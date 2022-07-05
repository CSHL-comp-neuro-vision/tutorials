% s = mkPlaid(stimSz, plaidDirection, gratingSf, gratingTf, plaidAngle,
%            plaidContrast)
%
% mkPlaid returns a drifting plaid stimulus
%
% Required arguments:
% stimSz            the size of the stimulus in [Y X T] coordinates
% plaidDirection    the direction of motion
% gratingSf         the spatial frequency of the plaid's grating components
% gratingTf         the temporal frequency of the plaid's grating components
%
% Optional arguments:
% plaidAngle        the angle between the grating components in radians.
%                   DEFAULT = (2/3)*pi (120 degrees).
% plaidContrast     the overall contrast of the plaid.

function s = mkSin(varargin)

% The following arguments are optional and by default are 'default'
plaidAngle = 'default';
plaidContrast = 'default';

% parse the varargin
                    stimSz = varargin{1};
                    plaidDirection = varargin{2};
                    gratingSf = varargin{3};
                    gratingTf = varargin{4};
if nargin >= 5      plaidAngle = varargin{5};       end
if nargin >= 6      plaidContrast = varargin{6};    end

% assign default values where appropriate
if strcmp(plaidAngle, 'default');       plaidAngle = (2/3)*pi;      end
if strcmp(plaidContrast, 'default');    plaidContrast = 1;          end

firstDirection = mod(plaidDirection + plaidAngle/2, 2*pi);
secondDirection = mod(plaidDirection - plaidAngle/2, 2*pi);
firstGrating =  mkSin(stimSz, firstDirection, gratingSf, gratingTf, ...
                      plaidContrast/2);
secondGrating = mkSin(stimSz, secondDirection, gratingSf, gratingTf, ...
                      plaidContrast/2);
s = firstGrating + secondGrating;
