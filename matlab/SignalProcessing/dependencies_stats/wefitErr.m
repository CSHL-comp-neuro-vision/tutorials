function err = wefitErr(X,levels,ncorrect,nerror,fixedBeta)
%
%AUTHOR:  Wandell
%DATE:  July 1, 1995
%PURPOSE:
%
%  Error function used in wefit computation.  This is the error to
% use to obtain a maximum likelihood fit of the Weibull through
% the data.  The derivation is in the appendix of a 1970s paper by Watson
% on some temporal measurements.
%
%ARGUMENTS:
%  X:  parameters passed by fmins
%  levels:  The contrast levels
%  ncorrect:  Number of correct responses
%  nerror:    Number of errors
%OPTIONAL:
%  fixedBeta: If beta is fixed in the search, pass it in.
%

if nargin == 5

  a = X(1);
  b = fixedBeta;

elseif nargin == 4
  a = X(1);
  b = X(2);

end

pr_correct = 1 - 0.5*exp(- ((levels/a).^b));

if (any(pr_correct==1))
 x = -Inf;
else
 x = sum(log10(1-pr_correct) .* nerror);
end

err = -( sum(log10(pr_correct) .* ncorrect) + x);


