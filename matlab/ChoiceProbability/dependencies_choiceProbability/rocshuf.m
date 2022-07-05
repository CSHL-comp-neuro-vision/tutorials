function [roc0, p] = rocshuf(x,y,nsims,dfname)
% [roc0, p] = rocshuf(x,y,nsims,dfname)
% 	Takes two distributions, x and y (column vectors) and returns the
% 	area under the ROC curve and the probability that this value could
% 	be obtained under the null hypothesis that x and y were draws from a
% 	common distribution, i.e., that the area under the ROC equals 0.5
% 	
% INPUT
% ~~~~~
% x,y   	responses.  Generally, x1 should have the larger mean (so
% 		the area comes out > 0.5
% nsims   	number of points in the bootstrap (2000 is a good idea)
% dfname 	string containing filename for matlab storage of variable A,
% 		the distribution of ROC areas under the null hypothesis. If
% 		a 4th argument is specified but is not a proper string, then
% 		the routine saves the distribution as rocshuftmp.mat
%
% Output
% ~~~~~~
% roc0	roc value using 100 criterion points.
% p	probability of getting a difference of roc vals under the
% 	hypothesis that x1 and x2 were drawn from a common parent
% 	distribution containing their union.  Under the null hypothesis, the
% 	size of x1 and x2 are fixed to the observed value.  That's why this
% 	is a shuffle.
%
% note: /home/monkeybiz/mike/mat/mfiles must be in the matlab path
%       relies on stats toolbox for randperm. see comment
% see also rocshuf2.m for comparing 2 roc areas.

% 5/2/97 update to matlab 5
% 1/24/98 no longer save the temp mat file unless there are 4 args

%%% test data
%$$$ x = randn(1,100);
%$$$ y = 1.5 * randn(100,1) + .2;
%$$$ nsims = 20

% t0 = clock;
z = [x(:);y(:)];			% parent distribution
b = logical([ones(size(x(:))); zeros(size(y(:)))]);
n0 = length(b);
roc0 = rocN(z(b),z(~b),100);
for i = 1:nsims
  I = randperm(n0);
  % [q I] = sort(rand(size(b)));
  B = b(I);
  a(i) = rocN(z(B),z(~B),100);
end
absa = abs(roc0 - .5);
p = sum(a >= absa+0.5)/nsims + sum(a <= 0.5-absa)/nsims;
if nargin > 3
  if ischar(dfname)
    s = sprintf('save %s.mat a',dfname)
    eval(s)
  else
    save /home/monkeybiz/mike/rocsim/soc/rocshuftmp.mat roc0 x y a
  end
end
% etime(clock,t0)				% print elapsed time

