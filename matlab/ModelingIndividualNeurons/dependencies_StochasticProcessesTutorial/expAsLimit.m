clear all

k = .1
m = [1:100]';
y = (1 - k./m).^m;

figure(1), clf
plot(m,y,'k-','LineWidth',3)
hold on
line(get(gca,'XLim'), exp(-k)*[1 1],'Color','r')
xlabel('m')
ylabel('y')


