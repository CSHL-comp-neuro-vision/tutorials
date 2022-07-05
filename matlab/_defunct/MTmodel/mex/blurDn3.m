% res = blurDn3(m, lev, filt)       Blur and downsample a 3D matrix
% 
% r = blurDn3(s); Blurs and downsamples s.
% r = blurDn3(s, 2)   does the same thing.
% r = blurDn3(s, 3)   does it twice.
% r = blurDn3(s, 2, filt)   uses the filter specified by filt to do the
% blurring. Filt is a row vector that specifies a separable 3D filter.

function res = blurDn3(m, lev, filt)

if exist('lev')~= 1
    lev = 2;
end
if exist('filt')~=1
    filt = [0.0884, 0.3536, 0.5303, 0.3536, 0.0884]';
end

if lev == 1
    res = m;
    return
end

fsz = length(filt);

res = validCorrDn3(m, reshape(filt, [1 1 fsz]));
res = validCorrDn3(res, reshape(filt, [1 fsz 1]));
res = validCorrDn3(res, reshape(filt, [fsz 1 1]));
res = res(1:2:end, 1:2:end, 1:2:end);

if lev > 2
    res = blurDn3(res, lev-1, filt);
end