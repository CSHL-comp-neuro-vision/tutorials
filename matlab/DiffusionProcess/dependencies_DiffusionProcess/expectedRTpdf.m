function RTpdfexp = expectedRTpdf(p,t,kRange)
% RTpdfexp = expectedRTpdf(p,[t],[kRange])

% Calculates the analytical solution for the distrubution RT's on correct
% trials for a diffusion model with the following parameters (as fields in
% the structure 'p'):

% p.a  upper bound (correct answer)
% p.b  lower bound (incorrect answer)
% p.s  standard deviation of drift (units/sec)
% p.u  drift rate (units/sec)

if ~exist('t')
    t = linspace(.1,30,101);
end

if ~exist('kRange')
    kRange = 5;
end

Pexp = expectedPC(p);

A = exp(-(t*p.u-2*p.a)*p.u/(2*p.s^2))./sqrt(2*pi*p.s^2*t.^3);
B = 0;

for k=-kRange:kRange
    B = B+ (p.a+2*k*(p.a+p.b))*exp(-(p.a+2*k*(p.a+p.b))^2./(2*t*p.s^2));
end

RTpdfexp = A.*B/Pexp;