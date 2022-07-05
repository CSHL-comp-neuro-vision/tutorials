function time_response=temp_imp_resp(n,k,t)
%time_response=temp_imp_resp(n,k,t)
%
%Produces a temporal impulse response function using the from from
%figure 1 in Adelson & Bergen (1985)
%
%It's pretty much a difference of Poisson functions with different
%time constants.

time_response=(k*t).^n .* exp(-k*t).*(1/factorial(n)-(k*t).^2/factorial(n+2));
