%StimulusSpace.m


n = 500;

[x,y] = meshgrid(linspace(0,1,n));
nt = 8;
dList = 30:30:360;
nd = length(dList);
posList = linspace(0,1,nd+3);
pos1 = posList(2)/2;
posList = posList(3:end);
g.r = .45;
sf = 5;
%g.sf = 30;  %cyclces/image
g.c = 1/6;

dListId = ceil(rand(1,6)*12);
%dListId = 
dListId = [1,1,1,1,1,1,1];
figure(1)
clf
hold on
for i=1:length(dList);
    nd = sum(dListId==i);
    if nd
       plot([0,nd*exp(sqrt(-1)*dList(i)*pi/180)],'ko-');
    end
end
axis equal
set(gca,'XLIm',[-3,3]);
set(gca,'YLim',[-3,3]);


figure(2)
clf
clear M
for i=1:nt
    img = zeros(n);
    g.ph = i*360/nt;

    g.xc = .5;
    g.yc = .5;
    for j=1:6
        g.d = dList(dListId(j));
        g.sf = sf;
        img = addGrating(img,g,x,y);
    end
    image( (img+1)*32);
    colormap(gray(64));
    axis equal
    axis off
    %showImage(img);
    set(gca,'YDir','normal');
    M(i) = getframe;
end

movie(M,20);
