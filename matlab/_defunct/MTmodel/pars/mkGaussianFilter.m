% res = mkGaussianFilter(sigma)         
% 
% Make a 1D gaussian filter with a given SIGMA. The size of the filter is 
% 3 * SIGMA

function res = mkGaussianFilter(sigma)

if sigma == -1
    res = 1;
else
    filterSize = ceil(3.*sigma);
    x = [-filterSize:filterSize];
    res = exp(-x.^2./(2.*sigma.^2));
    res = res./sum(res);
end