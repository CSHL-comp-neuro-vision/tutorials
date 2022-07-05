function drawEllipse(mu,cov,color)

%-------------------------------------------
%
% drawEllipse(mu,cov,color)
%
% draw an ellipse showing the 1 SD contour of a 
% two-dimensional gaussian
%
% required input:
% mu -- mean of the gaussian
% cov -- covariance matrix of the gaussian
%
% optional inputs:
% color -- line color (default = 'k')
%
% freeman, 06-23-2010
%-------------------------------------------

if ~exist('color','var');
    color = 'k';
end

theta = linspace(0,2*pi,100);
x = sin(theta);
y = cos(theta);

trans = [x', y']*sqrtm(cov);
trans = bsxfun(@plus,trans,mu);

hold on
plot(trans(:,1),trans(:,2),'color',color,'LineWidth',1);
hold off