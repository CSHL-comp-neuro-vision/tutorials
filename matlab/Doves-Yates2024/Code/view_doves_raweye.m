%Example to read raw eye movements directly from the .eye file in DOVES

% Author    : Umesh Rajashekar, Ian van der Linde
% Date      : 27.07.01, revised 15.05.04


clear;close all;

% Define Constants
IMAGES_DIR      = '../Images/';
DATA_DIR        = '../RawData/';


%% Select Input Files
[file_name DATA_DIR] = uigetfile([DATA_DIR '*.eye'],'Please select a Subject');


%% Process Fixation Information                                          %

image_count     = 0;
eye_fid = fopen([DATA_DIR file_name],'r');

while(~feof(eye_fid))

    [image_name,eye_x,eye_y,eye_pupil]  = read_doves_raweye(eye_fid);
    image_name                          = deblank(image_name);
    image_path                          = [IMAGES_DIR image_name];
    my_image                            = read_vanhat_foranalysis(image_path);

    if (length(eye_x>0))

        figure(1);clf;
        imagesc(my_image.^0.8);colormap gray;axis image; axis off; hold on;
        set(gca,'Position',[0 0 1 1]);

        h = plot(eye_x,eye_y,'Color',[.9 .9 .9],'Linewidth',2);axis ij;

        image_count = image_count + 1;
        display('Hit any key to continue');
        pause;


    end
end
fclose(eye_fid);


