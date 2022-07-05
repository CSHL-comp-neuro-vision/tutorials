function k = kurtMat(X)
%k = kurtMat(X)

X= X - repmat(mean(X),size(X,1),1);
num = mean(X.^4);
denom = mean(X.^2);

k = norm(num./denom.^2);


