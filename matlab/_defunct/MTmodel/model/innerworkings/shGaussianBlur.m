
function [res, ind] = shGaussianBlur(varargin);

pop = varargin{1};
popind = varargin{2};
nscales = size(popind, 1) - 1;
f = varargin{3};
fsz = length(f);
fx = reshape(f, [1 fsz 1]);
fy = reshape(f, [fsz 1 1]);
if nargin > 3
    ft = varargin{4};
    ftsz = length(ft);
    ft = reshape(ft, [1 1 ftsz]);
else
    ft = 1;
    ftsz = 1;
end

ind = shTrimIndices(popind, zeros([fsz fsz ftsz]));
res = zeros(ind(end, 1), size(pop, 2));
for s = 1:nscales
    for n = 1:size(pop, 2)
        tmp = shGetSubPop(pop, popind, n, s);
        
        tmp = validCorrDn3(tmp, fx);
        tmp = validCorrDn3(tmp, fy);
        if ft ~= 1
            tmp = validCorrDn3(tmp, ft);
        end
        res = shSetSubPop(res, ind, tmp, n, s);
    end
end
