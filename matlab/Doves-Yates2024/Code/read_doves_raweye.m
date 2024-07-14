function [image_name, eye_x,eye_y,eye_pupil] = read_doves_raweye(eye_fid)

%Function to read the raw eye movement data from .eye files in the
%DOVES database

% Author    : Umesh Rajashekar, Ian van der Linde
% Date      : 27.07.01, revised 15.05.04

%while(1)
while(~feof(eye_fid))
    temp = fgets(eye_fid);
    if (strfind(temp,'Eye Sample Data'))
        break;
    end
   
end

    if (feof(eye_fid))
        image_name = []; eye_x=[];eye_y=[];eye_pupil=[];
        return;
            
    end

    image_name = fgets(eye_fid);
    num_eye_points = fscanf(eye_fid,'%d\n',1);	%read y value of target
    eye_data = fscanf(eye_fid,'%f %f %f\n',num_eye_points*3);	%read the x y and pupil of eye data

    eye_x = eye_data(1:3:length(eye_data));
    eye_y = eye_data(2:3:length(eye_data));
    eye_pupil = eye_data(3:3:length(eye_data));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: write_eye_data                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteEyeData(eye_fid, current_image_name, temp_eye, current, total)
fprintf(eye_fid, '[Eye Sample Data (%d of %d)]\n', current, total);
fprintf(eye_fid, '%s\n', current_image_name);
fprintf(eye_fid, '%d\n', size(temp_eye, 2));
fprintf(eye_fid, '%f %f %f\n', temp_eye);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: write_task_data                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteTaskData(eye_fid, image_name, fixnumber, xpos, ypos, choice)
fprintf(eye_fid, '[Task Data]\n');
fprintf(eye_fid,'Image Used in Task: %s\n', image_name);
fprintf(eye_fid,'Fixation Number: %d\n', fixnumber);           % -1 if random
fprintf(eye_fid,'X Coordinate: %f\n', xpos);
fprintf(eye_fid,'Y Coordinate: %f\n', ypos);
fprintf(eye_fid,'User Response: %d\n', choice);
return