%StimulusSpace.m
%
%11/1/06 G.M. Boynton, Salk Institute

n = 500; %resolution of entire image (pixels)

[x,y] = meshgrid(linspace(0,1,n));
nt = 8;  %number of temporal frames/cycle
dList = 30:30:360; %direction list
nd = length(dList);

posList = linspace(0,1,nd+3);  
pos1 = .75*posList(2);
posList = posList(3:end);
g.r = .425/(nd+2);  %radius of each grating/plaid
sf = [70,70]; %spatial frequency of plaid components (cycles/image)
g.c = 1/2;  %contrast of each component.

figure(1)
clf
set(gca,'Position',[0,0,1,1]);

clear M
for i=1:nt
    img = zeros(n);
    g.ph = i*360/nt;
    for d1 = 1:nd
        for d2 = 1:nd
            g.xc = posList(d1);
            g.yc = posList(d2);

            g.d = dList(d1);
            g.sf = sf(1);
            img = addGrating(img,g,x,y);
            g.d = dList(d2);
            g.sf = sf(2);
            img = addGrating(img,g,x,y);
        end

        g.xc = posList(d1);
        g.yc = pos1;
        g.d = dList(d1);
        g.sf = sf(1);
        img = addGrating(img,g,x,y);

        g.xc = pos1;
        g.yc = posList(d1);
        g.d = dList(d1);
        g.sf = sf(2);
        img = addGrating(img,g,x,y);

    end



    image((img+1)*32);
    colormap(gray(64));
    axis equal
    axis off
    set(gca,'YDir','normal');
    M(i) = getframe;
end

movie2avi(M,'StimulusSpace');
%movie(M,40);
