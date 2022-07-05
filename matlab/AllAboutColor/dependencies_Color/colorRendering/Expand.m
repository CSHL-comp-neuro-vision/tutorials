function A=Expand(A,horizontalFactor,verticalFactor)
% A=Expand(A,horizontalFactor,[verticalFactor])
%
% Expands the matrix A by cell replication, and returns the result.
% If the vertical scale factor is omitted, it is assumed to be 
% the same as the horizontal. Note that the horizontal-before-vertical
% ordering of arguments is consistent with image processing, but contrary 
% to MATLAB's rows-before-columns convention.
%
% We use "Tony's Trick" to replicate a vector, as explained
% in MathWorks MATLAB Technote 1109, section 4.
%
% Also see ScaleRect.m
%
% Denis Pelli 5/27/96, 6/14/96, 7/6/96

if nargin<2 | nargin>3
	error('Usage: A=Expand(A,horizontalFactor,[verticalFactor])');
end
if nargin==2
	verticalFactor=horizontalFactor;
end
if round(verticalFactor)~=verticalFactor | verticalFactor<1 ...
	round(horizontalFactor)~=horizontalFactor | horizontalFactor<1
	error('Expand only supports positive integer factors.');
end
if isempty(A)
	error('Can''t expand an empty matrix');
end
[n,m]=size(A);
if horizontalFactor~=1
	A=reshape(A',1,n*m);
	A=A(ones(1,horizontalFactor),:);
	m=m*horizontalFactor;
	A=reshape(A,m,n)';
end
if verticalFactor~=1
	A=reshape(A,1,n*m);
	A=A(ones(1,verticalFactor),:);
	n=n*verticalFactor;
	A=reshape(A,n,m);
end
