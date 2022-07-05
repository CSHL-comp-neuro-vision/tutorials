function img = makeNeuralImage(p,pN)
%img = makeNeuralImage(p,pN)
%
%Creates a 'neural image' on the axes defined by p.x and p.theta
%based on the list of position and orientation centers and widths in
%'pN'.

if ~isfield(pN,'contrast')
    pN.contrast = ones(size(pN.x.center));
end

if ~isfield(pN,'range')
    pN.range = [0,1];
end

if ~isfield(pN,'method')
    pN.method = 'multiply';
end

x = p.x(:)';
theta = p.theta(:);

img= zeros(length(theta),length(x));
%loop through the list of stimulus components
for i=1:length(pN.x.center);
    if isfinite(pN.x.width(i))
        xG = pN.x.width(i)*sqrt(2*pi)*normpdf(x,pN.x.center(i),pN.x.width(i));
    else
        xG = ones(size(x));
    end
    if isfinite(pN.theta.width(i))
        thetaG =  pN.theta.width(i)*sqrt(2*pi)*normpdf(theta,pN.theta.center(i),pN.theta.width(i));
    else
        thetaG = ones(size(theta));
    end
    switch pN.method
        case 'multiply'
            subImg = thetaG*xG;
        case 'add'
            subImg = .5*repmat(thetaG,1,length(p.x))+.5*repmat(xG,length(p.theta),1);
        case 'or'
            thetaImg = repmat(thetaG,1,length(p.x));
            xImg = repmat(xG,length(p.theta),1);
            
            subImg = 1-(1-thetaImg).*(1-xImg);
    end
    img = img+pN.contrast(i)*subImg;
end

img = img*(pN.range(2)-pN.range(1))+pN.range(1);



