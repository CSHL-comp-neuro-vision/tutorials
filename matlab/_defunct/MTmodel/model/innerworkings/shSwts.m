% res = shSwts(dirs)        get weights for interpolating separable filters
% dirs is an mx3 matrix specifying the direction of m neurons in 3d
% fourier space. res is a mxn matrix, where n is the number of separable
% derivatives. The value of n depends on the order of the derivative. res
% contains a matrix of weights you can use to interpolate from the
% separable derivative responses you already have to the responses of other
% derivative filters (the ones pointing in the directions given by dirs).
%
% Example of use: say you have S (as returned by shModelV1Linear), a matrix
% containing the responses of separable derivative filters; say you also
% have dirs, as defined above. Then
% R = shSwts(dirs) * S;
% where R contains the responses of the neurons specified by dirs.

function [res, dirs] = shSwts(dirs)

% switch back to rectangular coordinates
dirs(:,2) = atan3(dirs(:,2), ones(size(dirs, 1), 1));
dirs = sphere2rec(dirs);

% normalize the direction vectors
d = sqrt(sum(dirs.^2, 2));
d = repmat(d, 1, size(dirs, 2));
dirs = dirs./d;

% precalculate factorials
for n = 1:3
    fac(n) = factorial(n);
end
fac = [1 fac];

% generate the weighting vectors
res = zeros(size(dirs,1), 10);
pt = 0;
for o3 = 0:3
    for o2 = 0:3-o3
        o1 = 3-o3-o2;
        pt = pt+1;
        const = fac(4)./(fac(o3+1)*fac(o2+1)*fac(o1+1));
        res(:,pt) = const .* dirs(:,1).^o1 .* dirs(:,2).^o2 .* dirs(:,3).^o3;
    end
end
