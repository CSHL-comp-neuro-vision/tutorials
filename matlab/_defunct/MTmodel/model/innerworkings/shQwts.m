% qWts = shqWts(dirs)      get the weights for interpolating squared directional derivative filters
% dirs is an mx3 matrix specifying the direction of m neurons in 3d
% fourier space; qWts is a matrix whose rows are the weights for
% interpolating from a population of qWtsponses that you already have.
%
% Example use: you have Q, N, or C. Assume for this example that you have
% Q. You also know the direction of the filters in Q, which we'll call
% v1dirs. You have dirs, which contains the directions of the
% neurons whose qWtsponses you want to get through interpolation. Let's call
% R the qWtsponses of these neurons.
% R = shqWts(dirs) * shqWts(v1dirs)^-1 * Q;
% That ought to do the trick.

function qWts = shQwts(dirs)

dirs(:,2) = atan3(dirs(:,2), ones(size(dirs, 1), 1));
dirs = sphere2rec(dirs);

d = sqrt(sum(dirs.^2, 2));
d = repmat(d, 1, size(dirs, 2));
dirs = dirs./d;

fac = zeros(1, 7);
fac(1:2) = 1;
for n = 2:6
    fac(n+1) = n*fac(n);
end

% generate the weighting vectors
qWts = zeros(size(dirs,1), 28);
pt = 0;
for o3 = 0:6
    for o2 = 0:6-o3
        o1 = 6-o3-o2;
        pt = pt+1;
        const = fac(7)./(fac(o3+1)*fac(o2+1)*fac(o1+1));
        qWts(:,pt) = const .* dirs(:,1).^o1 .* dirs(:,2).^o2 .* dirs(:,3).^o3;
    end
end