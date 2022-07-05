% S = mkFract(DIMS, VEL, FRACT_DIM, AMPL)     make a drifting fratal noise stimulus
%
% VEL is the velocity of the stimulus in [y, x] coordinates. It must be
% integers, because this hacky little program just uses circshift on
% fractal noise images generated using Eero's old code.
%
% DIMS is the size of the image
% FRACT_DIM is 1 by default.
% AMPL is the amplitude of the stimulus.

function res = mkFract(dims, vel, fract_dim, ampl)

if (exist('fract_dim') ~= 1)
    fract_dim = 1.0;
end

tdim = dims(3);
dims = dims(1:2);


res = randn(dims);
fres = fft2(res);

sz = size(res);
ctr = ceil((sz+1)./2);

shape = ifftshift(mkR(sz, -(2.5-fract_dim), ctr));
shape(1,1) = 1;  %%DC term

fres = shape .* fres;
fres = ifft2(fres);

if (max(max(abs(imag(fres)))) > 1e-10)
    error('Symmetry error in creating fractal');
else
    res = real(fres);
    res = res / sqrt(var2(res));
end

tmp = res;
res = zeros(dims(1), dims(2), tdim);
for t = 1:tdim
    tmp = circshift(tmp, [-vel(1) vel(2)]);
    res(:,:,t) = tmp;
end

res = res - min2(res);
res = ampl.*res./max2(res);
