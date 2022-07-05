% [pop, norm] = shModelMtNormalization(linpop, pars)  get squared and normalized linpop responses
% linpop is the response returned by shModelMtLinear; pars is the SH parameters
% structure given by, e.g. shPars. pop contains the normalized
% responses; norm contains the normalization signal.
% Extraneurons is a four-dimensional matrix containing the responses of
% neurons not in the the population

function [res, ind, pop] = shModelMtNormalization_Self(linpop, ind, pars, linres);

mtsigma = pars.mtsigma;

pop = linpop./(linpop + mtsigma);
res = linres./(linres + mtsigma);
