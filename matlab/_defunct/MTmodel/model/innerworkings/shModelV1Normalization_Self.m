% [res, ind, pop, resnorm, popnorm, Qres, Qpop] = shModelV1Normalization_Tuned1(Qpop, Qind, Qres, pars, extradirs)
% do "type 1" tuned normalization on V1 neurons
 

function [res, ind, pop, resnorm, popnorm] = shModelV1Normalization_Tuned1(Qpop, Qind, Qres, pars, extradirs)

v1sigma = pars.v1sigma;         % semi-saturation constant of normalization signal
v1dirs = pars.v1dirs;           % directions (in 3d Fourier space) of the v1 filters
scales = pars.scales;           % number of spatial scales
normstrength = pars.v1normstrength;

resnorm = normstrength.*Qres + v1sigma;
popnorm = normstrength.*Qpop + v1sigma;
res = Qres./(resnorm);
pop = Qpop./(popnorm);
