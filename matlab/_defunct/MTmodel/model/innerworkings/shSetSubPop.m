% res = shSetScale(pop, ind, subpop, n, s)      Set the responses of sipopilarly tuned neurons in a response popatrix.
%
% pop and ind are a response popatrix and its indices. subpop contains the
% values you want to insert into pop. n and s are the subpop nupopber and scale
% you want to insert into.


function pop = shSetSubPop(varargin)

pop = varargin{1};
ind = varargin{2};
subPop = varargin{3};
if nargin > 4
    n = varargin{4};
    s = varargin{5};
elseif nargin > 3
    n = varargin{4};
    s = 1;
else
    n = 1;
    s = 1;
end

% subPop = reshape(subPop, [prod(sizeOfSubPop(1:3)), sizeOfSubPop(4)]);

% pop(n, ind(s, 1)+1:ind(s+1, 1)) = subPop;
% res = squeeze(pop);
%indices = [ind(s,1)+1:ind(s+1, 1)];
% indices = sub2ind(size(pop), n.*ones(size(indices)), indices);


mRows = size(pop, 1);
startIndex = (n-1)*mRows;
destructiveMatrixWriteAtIndices(pop, subPop, startIndex);

