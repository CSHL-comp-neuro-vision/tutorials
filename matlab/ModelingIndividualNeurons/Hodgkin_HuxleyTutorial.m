% Hodgkin_Huxley.m
%
% Matlab implementation of the Hodgkin_Huxley model
% Uses the simplest (lowest order) approximation of the differential
% equation, so it's easy to read, but falls appart for large time steps 
%
% Written on 11/20/06 by G.M. Boynton at the Salk Institute, SNL-B
% Modified for CSHL 2014, 7/18/14

clear all; close all;

Iin_max = 10;   %driving (input) current (should be greater than about 6.5 mV to spike)
Iin_t = [30,150];

% Parameters for voltage clamp
Vclamp = NaN;   %voltage clamp desired level (NaN turns it off)
gClamp = 4;

%Timing parameters:
dt = .005;      %time steps (ms) (should be <=.01)
tMax = 200;     %max time (ms); (one spike lasts about 17ms)

%Equilibrium potentials
VNa = 115;      %sodium 
VK = -12;       %potassium 
Vl = 10.613;    %extracellular 

C=1;            %membrane capacitance
gL = 0.3;       %leak conductance
T=6.3;          %temperature (C)
k=3^(.1*T-.63); %constant (function of T)
t = 0:dt:tMax;  %time vector

%zero out vectors
V = zeros(1,length(t));  %membrane potential
n = zeros(1,length(t));  %'n' gate
m = zeros(1,length(t));  %'m' gate
h = zeros(1,length(t));  %'h' gate
gNa = zeros(1,length(t));  %Sodium conductance
gK  = zeros(1,length(t));  %Potassium conductance
Iin = zeros(1,length(t));  %Time-course of input current

%nonzero initial parameters (important!)
m(1)=0.0505;
n(1)=0.3835;
h(1)=0.4782;

%Generate time-course for input current. 
Iin(t>Iin_t(1) & t<Iin_t(2)) = Iin_max;

for i=1:length(t)-1    
    %voltage dependent alpha values for m,n and h
    am = (2.5-.1*V(i))/(exp(2.5-.1*V(i))-1);
    an = (1-.1*V(i))/(10*(exp(1-.1*V(i))-1));
    ah = .07*exp(-V(i)/20);

    %voltage dependent beta values for m,n and h
    bm = 4*exp(-V(i)/18);
    bn = 0.125*exp(-V(i)/80);
    bh = 1/(exp(3-.1*V(i))+1);

    %differential equations for m,n and h
    m(i+1) = m(i)+ dt*k*(am*(1-m(i))-bm*m(i));
    n(i+1) = n(i)+ dt*k*(an*(1-n(i))-bn*n(i));
    h(i+1) = h(i)+ dt*k*(ah*(1-h(i))-bh*h(i));

    %sodium and potassium conductances (functions of m,n and h)
    gNa(i+1) = 120*m(i+1)^3*h(i+1);
    gK(i+1) = 36*n(i+1)^4;

    %voltage clamp (if Vclamp isn't NaN)
    if ~isnan(Vclamp)
        %V(i) = Vclamp;
        Iin = Iin - gClamp*(V(i)-Vclamp);  %Alternate method: adjust current using feedback
    end
    %Differential equation for membrane potential
    V(i+1)=V(i)+dt/C*(Iin(i)-gNa(i+1)*(V(i)-VNa)-gK(i+1)*(V(i)-VK)-gL*(V(i)-Vl));
end

%The rest is stuff for plotting.  

fontSize = 10;
figure(1)  %Membrane potential over time
clf
plot(t,V,'LineWidth',2);
xlabel('Time (ms)','FontSize',fontSize);
ylabel('Voltage (mv)','FontSize',fontSize);
set(gca,'YLim',[-20,120]);
title(sprintf('%5.2g mV injected',max(Iin)),'FontSize',fontSize);
hold on
plot(t(Iin>0),-20*ones(1,sum(Iin>0)),'r-','LineWidth',4);
set(gcf,'PaperPosition',[1,1,11,2]);
set(gca,'FontSize',fontSize');

figure(2) %Na nd K conductance over time
clf
hold on
plot(t,gNa,'r-');
plot(t,gK,'g-');
legend({'gNa','gK'},'FontSize',fontSize);
xlabel('Time (ms)','FontSize',fontSize);
ylabel('Conductance','FontSize',fontSize);
set(gca,'FontSize',fontSize');

figure(3)  %m,n and h gates over time
clf
hold on
plot(t,m,'r-');
plot(t,n,'g-');
plot(t,h,'b-');
legend({'m','n','h'},'FontSize',fontSize);
xlabel('Time (ms)','FontSize',fontSize);
ylabel('p','FontSize',fontSize);
set(gca,'FontSize',fontSize');

figure(4)  %phase plots for m,n and h vs membrane potential
clf
subplot(1,3,1)
plot(m,V);
xlabel('m','FontSize',fontSize);
ylabel('V (mV)','FontSize',fontSize);
set(gca,'FontSize',fontSize');

subplot(1,3,2)
plot(n,V);
xlabel('n','FontSize',fontSize);
set(gca,'FontSize',fontSize');

subplot(1,3,3)
plot(h,V);
xlabel('h','FontSize',fontSize);
set(gca,'FontSize',fontSize');
