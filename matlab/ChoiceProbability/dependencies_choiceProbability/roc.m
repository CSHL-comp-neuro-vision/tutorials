function a = roc(x,y)
% ROC	computes area under ROC given distributions x and y
%	usage area = roc(x,y).  Defunct.  Use rocN.

a = rocN(x,y,100);
% $$$ 
% $$$ [m n] = size(x);
% $$$ x = reshape(x,1,m*n);
% $$$ [m n] = size(y);
% $$$ y = reshape(y,1,m*n);
% $$$ z= sort([x,y]);
% $$$ [m,n] = size(z);
% $$$ fa = zeros(1,n);	% allocate the vector
% $$$ hit = zeros(1,n);
% $$$ for i = n:-1:1
% $$$ 	fa(n-i+1) = sum(y > z(i));
% $$$ 	hit(n-i+1) = sum(x > z(i));
% $$$ end
% $$$ [m,ny] = size(y);
% $$$ fa = fa/ny;
% $$$ [m,nx] = size(x);
% $$$ hit = hit/nx;
% $$$ fa(1) = 0;
% $$$ hit(1) = 0;
% $$$ fa(n) = 1;
% $$$ hit(n) = 1;
% $$$ a = trapz(fa,hit);
% $$$ % uncomment next line if you want to see the plot
% $$$ % plot(fa,hit),axis('square'),xlabel('FA'),ylabel('Hit');


