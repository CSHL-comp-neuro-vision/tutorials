function y = sortXbyColI(x,i)
% y = sortXbyColI(x,i) sorts matrix, x, according to the order of column i.
% 	There are not provisions for ties.  The rows retain their integrity.
% 	The routine is necessary because applying matlab's sort to an array
% 	will sort each column independently.

[z I] = sort(x(:,i));
y = x(I,:);
