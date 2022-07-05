function ind = shTrimIndices(ind, filt)

fsz = [size(filt, 1), size(filt, 2), size(filt, 3)] - [1 1 1];

for i = 2:size(ind, 1)
    ind(i,2:4) = ind(i, 2:4) - fsz;
    ind(i, 1) = ind(i-1) + prod(ind(i,2:4));
end
