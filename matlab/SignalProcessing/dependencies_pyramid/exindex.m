function arr = exindex(arr, varargin)
%EXINDEX - extended array indexing
%   ARROUT = EXINDEX(ARRIN, S1, S2, ...) indexes a virtual array made by
%   extending ARRIN with zeros in all directions, using subscripts S1, S2
%   etc.
%
%   ARROUT = EXINDEX(ARRIN, S1, R1, S2, R2, ...) extends ARRIN using rule
%   R1 on its first dimension, R2 on its second dimension etc.
%
%   ARROUT = EXINDEX(ARRIN, S2, S2, ..., R) extends ARRIN using rule R on
%   every dimension.
%
%   Subscripts
%   ----------
%
%   Broadly, if V is the virtual extended array, ARROUT = V(S1, S2, ...)
%
%   The elements of S1, S2 etc must be integers. They need not be positive
%   and are not restricted in any way by the size of ARRIN. Logical
%   indexing and linear indexing are not supported.
%
%   There should be one subscript argument for each dimension of ARRIN,
%   unless ARRIN is a scalar or vector; in these cases one or two
%   subscripts may be used.
% 
%   The size of ARROUT is given by the normal Matlab rules for the result
%   of indexing into ARRIN: that is, size(ARROUT) = size(ARRIN, T1, T2,
%   ...) where T1 = ones(size(S1)) etc.
%
%   A subscript argument may be the string ':'. This behaves like a colon
%   in ordinary subscripting: a colon for the K'th subscript stands for
%   1:size(ARRIN, K). The 'end' keyword is not supported.
%   
%   Rules
%   -----
%
%   Each rule may be one of the following:
%
%   A number (a numerical scalar). The array is padded out with elements
%   equal to this value. A value of 0 is the default rule.
%   
%       If different constants are used on different dimensions, the
%       padding is done in the order of the subscripts. For example, a 2D
%       array is extended first in the row index direction and then in the
%       column index direction. For all other cases, the order in which
%       dimensions are extended has no effect.
%
%       There is a potential ambiguity if a vector or scalar is followed by
%       two arguments, and the second of these could be an subscript or a
%       rule. In this case, if the argument is intended to be a rule, the
%       value must be placed in a scalar cell; otherwise it will be treated
%       as an index.
%
%   'circular': ARRIN is extended with copies of itself; V is tiled with
%   ARRIN.
%
%   'symmetric': ARRIN is extended with copies of itself with reflection at
%   its boundaries; V is tiled with [ARRIN fliplr(ARRIN); flipud(ARRIN)
%   fliplr(flipud(ARRIN))].
%
%   'replicate': The array is padded out by copying its border elements; an
%   element of V is equal to the nearest element of ARRIN.
%
%   Examples
%   --------
%
%   Pad a 2D matrix with K extra rows and columns with reflection on both
%   axes:
%
%       b = exindex(a, 1-k:size(a,1)+k, 1-k:size(a,2)+k, 'symmetric');
%
%   Circularly shift a 2D matrix by R rows downward and C columns
%   rightwards:
%
%       b = exindex(a, 1-r:size(a,1)-r, 1-c:size(a,2)-c, 'circular');
%
%   Force a row or column vector to be 1024 elements long, trimming or
%   padding with zeros as necessary:
%
%       u = exindex(v, 1:1024);
%
%   or if the padding value was non-zero, say -1:
%
%       u = exindex(v, 1:1024, {-1});   % note constant in cell
%
%   If v was known to be a row vector, this could be written as:
%
%       urow = exindex(vrow, 1, 1:1024, -1);
%
%   See also: padarray, circshift

% Copyright David Young 2010

% Sort out arguments
[exindices, rules] = getinputs(arr, varargin{:});
consts = cellfun(@isnumeric, rules);  % Check for constants, as can be
constused = any(consts);              % more efficient if there are none

% Setup
nd = ndims(arr);
sz = size(arr);
if constused
    filled = cell(1, nd);
    c = zeros(1, nd);
end

% Main loop over subscript arguments, transforming them into valid
% subscripts into arr suing the rule for each dimension
for i = 1:nd
    ind = exindices{i};
    s = sz(i);
    if constused
        [exindices{i}, filled{i}, c(i)] = extend(ind, rules{i}, s);
    else % no need for information for doing constants
        exindices{i} = extend(ind, rules{i}, s);
    end
end

% Create the new array by indexing into arr. If there are no constants,
% this does the whole job
arr = arr(exindices{:});

% Fill areas that need constants
if constused
    % Get full range of output array indices
    ranges = arrayfun(@(x) {1:x}, size(arr));
    for i = nd:-1:1    % order matters
        if consts(i)
            ranges{i} = ~filled{i};     % don't overwrite original
            arr(ranges{:}) = c(i);      % fill with constant
            ranges{i} = filled{i};      % don't overwrite this
        end
    end
end

end

% -------------------------------------------------------------------------

function [exindices, rules] = getinputs(arr, varargin)
% Sort out and check arguments. Inputs are as given in the help comments
% for exindex. Outputs are cell arrays; each element of exindices is a
% set of integer extended indices which has been checked for validity; each
% element of rules is a rule which has not been checked for validity.

nd = ndims(arr);
nv = length(varargin);

if nd == 2 && min(size(arr)) == 1 && nv < 3
    % We have a vector or scalar ...
    if nv == 1  % ... with a single index and default rule
        if size(arr,1) == 1
            exindices = [{1} varargin];
        else
            exindices = [varargin {1}];
        end
        rules = {0 0};
    elseif nv == 2
        if isnumeric(varargin{2}) || strcmp(varargin{2}, ':')
            % ... with two indices and default rule
            exindices = varargin;
            rules = {0 0};
        else
            % ... with a single index and specified rule
            if size(arr,1) == 1
                exindices = [{1} varargin(1)];
            else
                exindices = [varargin(1) {1}];
            end
            rule = varargin{2};
            if iscell(rule) && isscalar(rule)
                % deal with special case of constant in cell
                rule = rule{1};
            end
            rules = {rule rule};
        end
    else    % in case nv = 0
        error('exindex:badnumargs', ...
            'Number of arguments inconsistent with array size');
    end
    
    % At least two indices
elseif nv == nd
    exindices = varargin;
    [rules{1:nd}] = deal(0);
elseif nv == nd+1;
    exindices = varargin(1:nd);
    [rules{1:nd}] = deal(varargin{end});
elseif nv == 2*nd
    exindices = varargin(1:2:end);
    rules = varargin(2:2:end);
else
    error('exindex:badnumargs', ...
        'Number of arguments inconsistent with array size');
end

% Expand any colons now to simplify checking.
% It's tempting to allow the 'end' keyword here: easy to substitute the
% size of the dimension. However, to be worthwhile it would be necessary to
% use evalin('caller',...) so that expressions using end could be given as
% in normal indexing. This would mean moving the code up to exindex itself,
% and evalin makes for inefficiency and fragility, so this hasn't been
% done.
colons = strcmp(exindices, ':');
if any(colons)  % saves a little time
    exindices(colons) = arrayfun(@(p) {1:size(arr,p)}, find(colons));
end

% Check the indices (rules are checked as required in extend)
checkindex = @(ind) validateattributes(ind, {'numeric'}, ...
    {'integer'}, 'exindex', 'index');
cellfun(checkindex, exindices);

end

% -------------------------------------------------------------------------

function [ind, filled, c] = extend(ind, rule, s)
% The core function: maps extended array subscripts into valid input array
% subscripts. 

if isnumeric(rule) && isscalar(rule) % pad with constant
    
    % The main messiness is due to constant padding. This can't be done
    % with indexing into the original array, but we want the indexing
    % structure to be preserved, so for now we index to element 1 on each
    % dimension, and record the indices of the regions that need to be
    % fixed.

    filled = ind >= 1 & ind <= s;
    ind(~filled) = 1;
    c = rule;
    
else    % pad with rule
    
    filled = [];  % never used
    c = 0;        % never used
    switch rule
        case 'replicate'
            ind = min( max(1,ind), s );
        case 'circular'
            ind = mod(ind-1, s) + 1;
        case 'symmetric'
            ind = mod(ind-1, 2*s) + 1;
            ott = ind > s;
            ind(ott) = 2*s + 1 - ind(ott);
        otherwise
            error('exindex:badopt', 'Unknown option');
    end
    
end

end


