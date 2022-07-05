function y = leakyIntegrator(s,k,t)

dt = t(2)-t(1);
nt = length(t);

y =zeros(1,nt);
for i=1:(nt-1)
    dy = s(i)-y(i)/k;
    y(i+1)=y(i)+dy*dt;
end