function S = powSpect(img, M, N)
% s=powSpect(img, M, N)
%	Calculates the power spectrum of the given img using a M x N
%	point 2D DFT.  If M and N are not specified the size of the image
%	is used instead.
%
% Rick Anthony
% 6/21/93

if (nargin < 3) [M N] = size(img); end

Y = fft2(img, M, N);
S = fftshift(Y.*conj(Y)/(M*N));
