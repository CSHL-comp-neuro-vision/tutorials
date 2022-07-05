function [alpha, beta, llik] = quickfit(data)
% QUICKFIT fits a weibull function to data using maximum likelihood
%	maximization under binomial assumptions.  It uses FITFUNW for
%	error calculation.  Data must be in 3 columns: x, %-correct, Nobs
%	usage: [alpha, beta] = quickfit(data)
%		alpha and beta are the threshold and slope parameters.

% update 4/1/95 to use optimization toolbox, contrained fit -- nope
%  old code is quickfit_no_opt.m
% 12/11/96  jdr & mns fixed guessing bug (we hope) 
% 06/25/06  GDF replaced fmins (now obsolete) with fminsearch.
% 06/17/10  MRC replaced table1 with interp1 
    


global Data;
Data = data;
% generate guess
q = ones(2,1);

% use linear interpolation to guess alpha
% next line takes x and %-cor columns, flips them and interpolates along
% x to find the value corresponding to .8 correct.  The interpolation
% requires monotonic %-cor column, so we sort the matrix 1st.
% Remember, it's just a guess.
a = find(data(:,2)>.7 & data(:,2)<.9);
% b = find(data(:,2)<.8 & data(:,2)>.8);
if isempty(a)
  q(1,1) = mean(data(:,1));
elseif min(data(:,2)) > 0.8
  q(1,1) = data(find(data(:,2)==min(data(:,2))),1);
elseif max(data(:,2)) >= 0.8
  % this bombs if there's nothing greater than 0.8 in the 2nd column
  tmpdata = data + [zeros(size(data(:,1))) rand(size(data(:,2)))/100 ...
	zeros(size(data(:,3)))];
  tmp = sortXbyColI(fliplr(tmpdata(:,1:2)),1);
  t = tmp(:,1); 
  p = tmp(:,2);
  q(1,1) = interp1(t,p,.8);
else
  q(1,1) = max(data(:,2));
end

trace = 0;
tol = .0001;
% quick = leastsq('qe2',q,[trace tol]);
options = optimset;
% foo = fitfunw([trace,tol]);
quick = fminsearch('fitfunw',q,options);
if quick(1) <= 0
  % IT'S A bad fit.  Try another method.
  q(1) = .1;
  q(2) = .5;
  quick = fmins('fitfunw',q);
  %$$$      vlb = [.000001 .001];
  %$$$      vub = [ ];
  %$$$      options(1) = 1;
  %$$$      options(12) = 0;
  %$$$      options(13) = 0;
  %$$$      [quick,g] = constr('quick_err', q, options, vlb, vub);
  
  %$$$   OPTIONS = 0;
  %$$$   OPTIONS(2)=1e-5;
  %$$$   OPTIIONS(3) = 1e-20;
  %$$$   OPTIONS(5) = 1;
  %$$$   OPTIONS(13) = 0;
  %$$$   OPTIONS(12) = 2;
  %$$$   [quick,g] = constr('quick_err', q, OPTIONS);
end
% quick = fmins('qe3',q,[trace tol]);
% quick = fmins('quick_err',q,[trace tol]);
% vlb = [.000001 .001];
% vub = [ ];
%$$$ options(1) = 1;
%$$$ options(12) = 2;
%$$$ options(13) = 0;
%$$$ % [quick,g] = constr('quick_err', q, options, vlb, vub);
%$$$ [quick,g] = constr('quick_err', q, options);
% quick
alpha = quick(1,1);
beta = quick(2,1);
llik = fitfunw(q);
 

