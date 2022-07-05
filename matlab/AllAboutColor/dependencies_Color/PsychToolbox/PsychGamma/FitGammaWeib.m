function [fit_out,x,err] = FitGammaWeib(values_in,measurements,values_out,x0)% [fit_out,x,err] = FitGammaWeib(values_in,measurements,values_out,x0)%% Fits a Weibull function to the passed gamma data.%% 10/3/93    dhb  Created from jms code to fit psychometric functions% Compute fitoptions = foptions;options(1) = 1;%options(14) = 600;vlb = [1e-5 1e-5]';vub = [1e5  10.0]';x = constr('FitGammaWeibFun',x0,options,vlb,vub,[],values_in,measurements);% Now compute fit values and error to data for returnfit_out = ComputeGammaWeib(x,values_out);fit_in = ComputeGammaWeib(x,values_in);err = ComputeRMSE(fit_in,measurements);