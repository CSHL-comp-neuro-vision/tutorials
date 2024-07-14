function [display_image] = read_vanhat_foranalysis(filename)

%Function to load in an image for the DOVES experiment

% Author    : Umesh Rajashekar, Ian van der Linde

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Define Constants                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

w = 1024; h = 768;                                              % W x H of Images(cropped van Hateren) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read van Hateren File                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1              = fopen(filename,'rb','ieee-be');
buf             = fread(f1,[w,h],'uint16');                     % Read 16bpp image data
display_image   = buf';                                         % Note: van Hateren images are flipped
display_image   = display_image+1;
fclose(f1);
display_image = display_image/max(display_image(:))*255;

return

% avg_grey = mean(display_image(:));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  display_image = floor(display_image/max(display_image(:)) * 6279);
%     
%     U               = unique(display_image);                        % Get the unique gray scales
%     H               = zeros(1,max(U)); H(U) = 1; H = cumsum(H);     % Mapping grayscales to numbers from 0-256(length of U)
%     new_map         = repmat(floor( length(U)* U/max(U)),1,3);    % Create a new lookup table
%     display_image   = H(display_image);                             % Create an indexed image
%     
