function result = isSingular(matrix)
% ISSINGULAR: Returns 1 if matrix is singular, 0 otherwise.
% 
% result = isSingular(matrix)
% 

if ((1/cond(matrix))<1.0e-6)
  result = 1;
else
  result = 0;
end  
  
