function [F,t] = real2complex(Y)
%[F,t] = real2complex(Y)
%
%Returns the complex-valued vector F in the convention of the output of the
%fft, given the real-valued vectors containing phases and amplitudes in the
%structure Y.  This is the inverse function of 'complex2real'
%
%Inputs:    Structure Y with fields:
%   dc       mean value of y
%   amp      vector of amplitudes (length ceil(length(t)/2))
%   ph       vector of phases (in degrees, cosine phase)
%   nt       length of t (needed for myifft)

%Outputs:
%   Y        complex-valued vector in the convention of the output of fft
%   t        time vector of size y (default is 1:length(y));
%
%SEE ALSO    complex2real fft ifft

%4/15/09     Written by G.M. Boynton at the University of Washington



F = zeros(1,Y.nt);

%dc
F(1) = Y.dc*Y.nt;

amp = Y.nt*Y.amp/2;
ph = -pi*Y.ph/180;

z = amp.*exp(sqrt(-1)*ph);

F(2:(length(Y.freq)+1)) = z;

F(end:-1:(end-length(z)+1)) = conj(z);

maxt = Y.freq(1);
dt = maxt/Y.nt;
t = 0:dt:(maxt-dt);