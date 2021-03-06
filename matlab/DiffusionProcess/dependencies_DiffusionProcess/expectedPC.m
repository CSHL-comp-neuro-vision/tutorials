function Pexp = expectedPC(p)
% Pexp = expectedPC(p)
%
% Calculates expected probability correct for a diffusion model with the
% following parameters (as fields in the structure 'p'):
% p.a  upper bound (correct answer)
% p.b  lower bound (incorrect answer)
% p.s  standard deviation of drift (units/sec)
% p.u  drift rate (units/sec)
%
% It's a one line function with the calculation:
% Pexp = (exp(2*p.u*p.b/p.s^2)-1)/(exp(2*p.u*p.b/p.s^2)-exp(-2*p.u*p.a/p.s^2));
Pexp = (exp(2*p.u*p.b/p.s^2)-1)/(exp(2*p.u*p.b/p.s^2)-exp(-2*p.u*p.a/p.s^2));
