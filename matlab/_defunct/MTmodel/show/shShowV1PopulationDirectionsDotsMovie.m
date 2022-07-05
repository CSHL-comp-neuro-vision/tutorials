% dirsMovie = shShowV1PopulationDirectionsDotsMovie(dirsOverTime, axisView)		Show a movie of the V1 neuron directions being chosen.
%
% The output dirsMovie is a MATLAB movie you can watch with the MOVIE
% function. Warning: this can take a really long time.
%
% Required arguments:
% dirsOverTime          the second output from shParsV1PopulationDirections
% 
% Optional arguments:
% axisView              a valid argument for a call to the VIEW function.
%                       You will view the sphere from this perspective.
%                       DEFAULT = 3
% Output:
% dirsMovie             a standard MATLAB movie you can view with MOVIE
%
% Example of use:
% [v1Neurons, dirsOverTime] = shParsV1PopulationDirections(28, 300);
% dirsMovie = shShowV1PopulationDirectionsDotsMovie(dirsOverTime);
% movie(dirsMovie);

function dirsMovie = shShowV1PopulationDirectionsDotsMovie(varargin);

axisView = 'default';

dirsOverTime = varargin{1};
if nargin >= 2;     axisView = varargin{2};         end

if strcmp(axisView, 'default');     axisView = 3;       end;





nsteps = size(dirsOverTime, 3);
np = size(dirsOverTime, 2);
xmins = min(min(dirsOverTime,[],2),[],3);
xmaxs = max(max(dirsOverTime,[],2),[],3);

h = gcf;

[x,y,z] = sphere(50);
for n = 1:nsteps
    
    
    clf
    hold on
    
    for p = 1:np/2
        plot3(dirsOverTime(1,p,n), dirsOverTime(2,p,n), dirsOverTime(3,p,n), 'ro', 'markerfacecolor', 'r');
    end
    for p = np/2+1:np
        plot3(dirsOverTime(1,p,n), dirsOverTime(2,p,n), dirsOverTime(3,p,n), 'ro', 'markerfacecolor', 'r');
    end

    
    axis([-1 1 -1 1 -1 1]);
    surfl(x,y,z);    
    colormap(gray);
    shading interp;

    %     alpha(.6);
    
    axis off
    axis equal
    axis vis3d
    view(axisView);     % set the angle you want here
    hold off;
   
    
    
    dirsMovie(n) = getframe;

end
