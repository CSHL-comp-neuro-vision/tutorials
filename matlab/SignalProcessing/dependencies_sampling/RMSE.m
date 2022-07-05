function out = RMSE( A, B)
% RMSE: root mean-squared-error between two matrices/images.
% 
% result=RMSE(A,B)
%
out = sqrt(mean(mean(abs(A-B).^2)));
