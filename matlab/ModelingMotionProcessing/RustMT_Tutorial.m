%% RustMT.m

addpath('dependencies_RustMTModel');
% Implementation of model from:

% Rust, N.C., Mante, V., Simoncelli, E.P., & Movshon, J.A. (2006). How MT
% cells analyze the motion of visual patterns. Nat Neurosci, 9 (11),
% 1421-1431.
%
% Requires folder 'RustMTModel' to be added to the path
% June 2010, gmb wrote it.

clear all
%% define model parameters
p.pList = 30:30:360;  %list of each V1 neuron's preferred directions
s.dList = 30:30:360;  %list of possible stimulus directions

p.b = pi;           %tuning width of V1 neurons
p.s1 = .1;          %'untuned' normalization factor
p.s2 = .1;          %'self' normalization factor

p.A = 1;            %MT scaling nonlinearity
p.B = 2;            %MT exponent nonlinearity

% make up some weights of V1 influence on MT response:

%use this line for pattern cell:
p.w = exp(-(linspace(-1,1,length(p.pList)).^2)/.4)-.3;

%use this line for a component cell:
%p.w = exp(-(linspace(-1,1,length(p.pList)).^2)/.01); 

%% Plot the linear weights of V1 neurons
figure(1)
clf
plot(p.pList,p.w);
xlabel('Preferred orientation (deg)');
ylabel('Weight');

%% Now we will make the 144 stimuli used to generate the matrix, M, which is the 
%response to each pairwise combination of gratings that make all possible 
%plaids (as in their figure 4).
count = 0;
s.c = zeros(length(s.dList)^2,length(s.dList));
for i=1:length(s.dList)
    for j=1:length(s.dList);
        count = count+1;
        s.c(count,i) = s.c(count,i)+1/6;
        s.c(count,j) = s.c(count,j)+1/6;
    end
end

% get the model responses to the 144 stimuli
M = predRustMT(p,s);
%reshape it into a square
M = reshape(M,length(s.dList),length(s.dList));

%% Display the matrix like Rust et al. do in figure 4
figure(2)
clf
imagesc(s.dList,s.dList,M);
colormap(gray)
hold on
contour(s.dList,s.dList,M,8,'-','LineWidth',2,'Color',[1,.5,.5]);
set(gca,'YDir','normal');
axis equal
axis tight
xlabel('Direction 1 (deg)');
ylabel('Direction 2 (deg)');
set(gca,'XTick',[0,180,360]);
set(gca,'YTick',[0,180,360]);


%% Response to gratings vs. 120 deg plaids
s.c= zeros(length(s.dList),length(s.dList));
for i=1:length(s.dList)
    s.c(i,i) = 1/6;
end

M = predRustMT(p,s);
figure(3)
clf
plot(s.dList,M,'r-')
hold on

ni = round(60/(s.dList(2)-s.dList(1)));

s.c= zeros(length(s.dList),length(s.dList));
for i=1:length(s.dList)
    s.c(i,mod(i-ni-1,length(s.dList))+1) = 1/6;
    s.c(i,mod(i+ni-1,length(s.dList))+1) = 1/6;
end

M = predRustMT(p,s);
plot(s.dList,M,'b-') 
ylim = get(gca,'YLim');
set(gca,'XTick',0:45:360);
set(gca,'XLim',[0,360]);
xlabel('Orientation');
ylabel('Response');

