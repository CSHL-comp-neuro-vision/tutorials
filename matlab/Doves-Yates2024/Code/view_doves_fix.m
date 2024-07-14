%Example to read and display raw eye movements and fixations from the .mat 
%files in Fixations folder in DOVES

% Author    : Umesh Rajashekar, Ian van der Linde
% Date      : 27.07.01, revised 15.05.04

clear; close all;

IMAGES_DIR            = '../Images/';
FIXATIONS_DIR = '../Fixations/';

%% Get a list of all images
temp_names                  = getfilenames(IMAGES_DIR,'.iml');	
for file_num = 1:size(temp_names,1)
    temp                    = temp_names(file_num,:);
    temp                    = deblank(temp);
    img_filenames_list{file_num}  = temp;
end

img_filenames_list = sort(img_filenames_list);

%% 
for i = 1:length(img_filenames_list)

    image_name = img_filenames_list{i}; %stored as imk00031.iml
    
    %Load  the image
    my_image = read_vanhat_foranalysis([IMAGES_DIR image_name]);

    %Load the fixations from all subjects for that images
    load ([FIXATIONS_DIR image_name '.mat']); %loads subj_names_list, fix_data, eye_data
    
    figure(1); clf;
    imagesc(my_image.^0.3);colormap gray;axis image; axis off; hold on;
    set(gca,'Position',[0 0 1 1]);

    %Display the fixations for each subject
    for iSubj = 1:length(subj_names_list) % 
        
        %Raw eye movement data for a subject
        eyeData = eye_data{iSubj}; 
        eyeX = eyeData(1,:); eyeY = eyeData(2,:);

        %The computed fixations for the subject
        fixData = fix_data{iSubj};
        fixX = fixData(1,:); fixY=fixData(2,:);fixDur = fixData(3,:);
        
        h1 = plot(eyeX, eyeY, 'Color', [.9 .9 .9], 'LineWidth',2); axis ij; 
        h2 = plot(fixX,fixY,'ko','MarkerFaceColor','g','MarkerSize',8); axis ij;
        h3 = plot(fixX(1),fixY(1),'sr','MarkerFaceColor','r','MarkerSize',12); axis ij;
        
        display('Hit any key to continue');
        pause; 
        delete(h1); delete(h2); delete(h3);
    end
    
end

%% 