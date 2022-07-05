function  img = addGrating(img,g,x,y)

id  = (x-g.xc).^2 + (y-g.yc).^2 <= g.r.^2;
ramp = x.*cos(g.d*pi/180)+y.*sin(g.d*pi/180);
img(id) = img(id) + g.c*cos(2*pi*g.sf*ramp(id)-g.ph*pi/180);

%img(id) = img(id)+1;
