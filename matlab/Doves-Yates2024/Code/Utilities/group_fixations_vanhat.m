%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Title     : Group Fixations & Times on van Hateren Images              %
% Author    : Umesh Rajashekar, Ian van der Linde                       %
% Date      : 27.05.04                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input  - The .eye files from the vanahat experiment
%Output - Data structure stored in GROUPED_EYE_DATA
%       - Each .mat file is stored as image_name.iml.mat.
%       - Loading this .mat files gives two cell arrays -
%       - 1. subj_names_list - This Cell array has subject initials E.g. - {'IVDL2';'UR2'; ...}
%       - 2. fix_data         - This Cell array has [fix_x;fix_y;fix_duration]
%       - 3. eye_data        - This Cell array has [eye_x; eye_y]

%8 Jun 2004 - Changed file to do auto update. If you want to force a
%subject's update, select it manually

clear;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Define Constants                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

EYE_DATA_DIR            = '~/urTmp/Doves/RawData/';
GROUPED_EYE_DIR         = '~/urTmp/Doves/Fixations/';

IMAGE_WIDTH         = 1024;         % Should read from file!!!
IMAGE_HEIGHT        = 768;
DIST2SCREEN         = 134;          
XDPI                = 64;           
YDPI                = 64;
SAMPLE_RATE         = 200;
NUM_FIX             = SAMPLE_RATE/10;

% Select All the Eye Files For Grouping
%======================================
file_count = 0;
temp_names = getfilenames(EYE_DATA_DIR,'eye');
for file_num = 1:size(temp_names,1)
    temp = temp_names(file_num,:);
    temp = deblank(temp);
    eye_filename_list{file_num}    = temp;
end
  
eye_filename_list = sort(eye_filename_list); % Sorting the observers in alph. order - no reason
gp_filenames_list=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the Fixations from the eye files and stored them in the grouped data structure%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for file_cnt = 1:length(eye_filename_list)
    fprintf(['Working on ' num2str(file_cnt) ' of ' num2str(length(eye_filename_list)) '\n']);
    eye_filename  = eye_filename_list{file_cnt}
    subj_name      = eye_filename;
    eye_fid = fopen([ EYE_DATA_DIR eye_filename],'r');
    image_count = 0; % Update for the impatient user
    while(~feof(eye_fid))
        
        %Read in the eye scan path
        [image_name,eye_x,eye_y,eye_pupil]  = read_vanhat_fixations(eye_fid);
        image_name = deblank(image_name);
        if (isempty(image_name))break; end; %read_vanhat_fixations returns [] if EOF
       
%       image_count = image_count + 1        %Update for the impatient user
        
        if (length(eye_x>0))
            
            %convert eyedata into degrees
            eye_x = eye_x'; eye_y = eye_y';
            
            %remove pixels values that lie outside the area of interest
            [eye_x eye_y] = eye_cleanup(eye_x,eye_y,eye_pupil,IMAGE_WIDTH,IMAGE_HEIGHT);            
            eye_x_deg     = pix2deg(eye_x, DIST2SCREEN, XDPI);eye_y_deg = pix2deg(eye_y, DIST2SCREEN, YDPI);            
            [fix_x_deg fix_y_deg fix_start_index fix_duration saccade_begin]    = doves_fix_dll(eye_x_deg,eye_y_deg,NUM_FIX);
            
            fix_x = deg2pix(fix_x_deg, DIST2SCREEN, XDPI); fix_y = deg2pix(fix_y_deg, DIST2SCREEN, YDPI);   
            
            if (length(fix_x>0))
                
                %Using image_name, load the grouped eye data for the image (if it exists)
                gp_filename = [image_name '.mat'];
                if (sum(strcmp(gp_filenames_list,gp_filename)))
                    load ([GROUPED_EYE_DIR gp_filename]); % Data loaded in variable - subj_names_list; fix_data;eye_data
                    subj_index = length(subj_names_list)+1;
                else % could not find image_name.mat grouped file - new image
                    subj_index = 1;
                    subj_names_list = {};
                    fix_data={};
                    eye_data={};
                end
                
                subj_names_list{subj_index} = subj_name;
                fix_data{subj_index}= [fix_x ;fix_y;fix_duration];
                eye_data{subj_index} = [eye_x;eye_y];
                
                %Sort the data in alphabetical order
                %[subj_names_list, sort_order]= sort(subj_names_list);
                %fix_data = fix_data(sort_order);
                
                save([GROUPED_EYE_DIR gp_filename],'subj_names_list','fix_data','eye_data'); %Save the variables in image_name.mat file
            end                        
        end        
        
    end	        % while(~feof(eye_fid))
    fclose(eye_fid);    
   
    %Update the list of .mat found - this is needed for the first time
    %update of all files
    temp_names                  = getfilenames(GROUPED_EYE_DIR,'.mat');	
    for file_num = 1:size(temp_names,1)
        temp     = deblank(temp_names(file_num,:));
        gp_filenames_list{file_num}  = temp;
    end
    
end	% for file_cnt = 1:length(eye_filename_list)