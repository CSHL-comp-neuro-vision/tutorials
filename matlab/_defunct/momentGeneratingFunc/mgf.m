function y = mgf(theta,x,pdf)
% y = mgf(theta,x,pdf) returns the moment generating function evaluated at
% theta. User supplies the pdf, such that plot(x,pdf) would represent the pdf.
% User should ensure that the integral of the pdf over the range of x is 1.

% 7/10/03 mns wrote it

theta = theta(:);
x = x(:);
pdf = pdf(:);
[X H] = ndgrid(x, theta);
[P H] = ndgrid(pdf,theta);
y = trapz(x, exp(H.*X)  .* P);
y = y(:);

    

