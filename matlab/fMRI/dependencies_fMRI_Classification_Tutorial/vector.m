function vec = vector(mtx)

%-------------------------------------------
%
% vector(mtx)
%
% vectorize the elements of mtx
% same as mtx(:)
%
% required input:
% mtx -- matrix, any dimensionality
%
% freeman, 09-25-2010
%-------------------------------------------

vec = mtx(:);
