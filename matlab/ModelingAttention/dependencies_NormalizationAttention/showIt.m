function showIt(p,img,range)

if ~exist('range','var')
    range = [min(img(:)),max(img(:))];
end

clf
image(p.x,p.theta,255*(img-range(1))/(range(2)-range(1)));
set(gca,'YTick',[-180:90:180]);
set(gca,'XTick',[-200:100:200]);
colormap(gray(255));
xlabel('Position');
ylabel('Orientation');
axis equal
axis tight