function [alpha, beta] =  wefit(levels,ncorrect,nerror,fixedBeta)
%
%  [alpha, beta] =  wefit(levels,ncorrect,nerror,fixedBeta)
%
%AUTHOR:  Wandell
%DATE:    June 29, 1995
%PURPOSE:
%  Find the best fitting Weibull parameters to a data set.  The Weibull
%  is fit assuming a 2IFC model, that is
%
%          pCorrect = 1 - 0.5 exp( - (level/alpha) .^ beta )
%
% The fit is based on the method described in Watson's paper from the 70s,
% using a maximum likelihood fit.
%
%ARGUMENTS:
%  levels:  The contrast levels of the stimuli
%  ncorrect:  The number of correct responses at each level
%  nerror:    The number of errors at each level
%
%OPTIONAL
%  fixedBeta:  To fit the data with beta fixed, pass along the
%              value of beta
%

alpha = mean(levels);

if nargin == 4

  X = alpha;

  %One dimensional search
  alpha = fmins('wefitErr',X,foptions,[],levels,ncorrect,nerror,fixedBeta);

else

  %Initialize both arguments
  beta = 2.0;
  X = [alpha beta];

  %Search
  X = fmins('wefitErr',X,foptions,[],levels,ncorrect,nerror);

  %Return
  alpha = X(1);
  beta = X(2);

end


