function sd=QuestSd(q)% sd=QuestSd(q)% Get the sd of the threshold distribution.% See QuestDemo,QuestCreate, QuestUpdate, QuestQuantile, QuestMean, QuestMode,% QuestSd, and QuestSimulate.%% Denis Pelli, 6/8/96% 3/1/97 dgp updated to use MATLAB 5 structs.% 4/12/99 dgp dropped support for MATLAB 4.% Copyright (c) 1996-1999 Denis Pelliif nargin~=1	error('Usage: sd=QuestSd(q)')endif length(q)>1	sd=zeros(size(q));	for i=1:length(q(:))		sd(i)=QuestSd(q(i));	end	returnendp=sum(q.pdf);sd=sqrt(sum(q.pdf.*q.x.^2)/p-(sum(q.pdf.*q.x)/p).^2);