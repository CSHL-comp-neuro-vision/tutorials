function logy2raw(base,precision)
%logy2raw([base],[precision])
%
%Converts X-axis labels from log to raw values.
%base:    	base of log transform; default base is e.
%precision:	number of decimal places
%
%Example:
% x=linspace(-3,0,11);
% plot(log(x),log(x.^2));
% logx2raw
% logy2raw

%SEE ALSO;   Logy2raw
%11/17/96	gmb	Wrote it.
%6/6/96	        gmb added precision argument
%01/30/02       gmb updated it to use cell arrays, and to use original
%                xtick values instead of converting labels.  This way,
%		multiple calls to this function doesn't keep converting
%		the axis.

if nargin==0
    base=exp(1);
end

origYTick = get(gca,'YTick');
newYTick = base.^(origYTick);

for i=1:length(newYTick)
    if exist('precision','var')
        newYLabel{i}= num2str(newYTick(i),sprintf('%%5.%df',precision));
    else
        newYLabel{i} = num2str(newYTick(i));
    end
end
set(gca,'YTickLabel',newYLabel);



