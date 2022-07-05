function err = fitfunw(q)
%FITFUNW Used by QUICKFIT.
%	FITFUNW(lambda) returns the error between the data and the
%	values computed by the current function of weibul params.
%	FITFUNW assumes a function of the form
%
%	  y = 1 - .5 * exp( -(x/q(1))^q(2) )
%
%	thus q(1) is alpha, and q(2) is beta.
%	The data is in columns such that Data(:,1) is abscissa 
%	Data(:,2) is observed percent correct (0..1)
%	Data(:,3) is number of observations.
%	The value of err is the -log likelihood of obtaining Data
%	given the parameters q.
%
global Data Plothandle

TINY = eps;
x = Data(:,1);
y = Data(:,2);
n = Data(:,3);
z = 1 - .5 * exp( -(x/q(1)).^q(2) );
z = z - eps.*(z == 1);
z = z + eps.*(z == 0);

llik = n .* y .* log(z) +  n .* (1-y) .* log(1-z);
% for moment, just return the square error
% err = norm(z-y);
err = -sum(llik);

