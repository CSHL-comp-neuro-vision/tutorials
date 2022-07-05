function cost = logLikelihoodROC(p,isSignal,response)

critRange = [-inf,p.criterion,inf];

for i=1:length(critRange)-1;
    pResp(i,1) = normcdf(critRange(i+1),p.noiseMean,p.sd)-...
        normcdf(critRange(i),p.noiseMean,p.sd);
    pResp(i,2) = normcdf(critRange(i+1),p.signalMean,p.sd)-...
        normcdf(critRange(i),p.signalMean,p.sd);
end

cost = 0;

for i=1:length(isSignal);
    cost = cost-log(pResp(response(i),isSignal(i)+1));
end
