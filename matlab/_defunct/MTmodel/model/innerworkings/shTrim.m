% [res, ind] = shTrim(L, Lind, trimmer)      Trim the edges off a matrix of model responses
%
% L is the response matrix you want trimmed, Lind is its index matrix, and
% trimmer is a 3-vector contianing the number of pixels you want trimmed
% off on both sides of each dimension. So you'll actually trim away twice
% the number of pixels in trimmer.

function [pop, ind] = shTrim(bigPop, bigInd, trimmer);

x = trimmer(1);
y = trimmer(2);
t = trimmer(3);

% Trim down the indices
ind = bigInd;
for i = 2:size(bigInd, 1)
    ind(i,2) = ind(i,2) - 2*x;
    ind(i,3) = ind(i,3) - 2*y;
    ind(i,4) = ind(i,4) - 2*t;

    ind(i,1) = ind(i-1, 1) + prod(ind(i,(2:4)));
end
pop = zeros(ind(end, 1),size(bigPop, 2));

% Trim down the population response
w = [];
for i = 1:size(bigInd, 1)-1
    tmp = ones(bigInd(i+1, 2:end));
    tmp(x+1:end-x, y+1:end-y, t+1:end-t) = 0; 
    tmp = tmp(:);
    w = [w, find(tmp)' + bigInd(i, 1)];
end
bigPop(w, :) = [];
pop = bigPop;
