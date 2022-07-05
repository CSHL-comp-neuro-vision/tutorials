function plotResp(t,s,y,col)

if ~exist('col','var')
    col = {'r','b'};
end

subplot(2,1,1)
plot(t,s,'-','LineWidth',2,'Color',col{1});
xlabel('Time (s)');
title('Input');
set(gca,'XLim',[min(t)-.05,max(t)+.05]);

subplot(2,1,2)
plot(t,y,'-','LineWidth',2,'Color',col{2});
xlabel('Time (s)')
title('Output');
set(gca,'XLim',[min(t)-.05,max(t)+.05]);


