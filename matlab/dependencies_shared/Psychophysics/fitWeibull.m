function logLikelihood = fitWeibull(p,results)


pred = Weibull(p,results.intensity);

pred = pred*.99;

logLikelihood = -sum(log(pred).*results.response) - sum(log(1-pred).*(1-results.response));

