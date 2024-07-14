function [orig_calib_locs, record_calib_locs] = read_vanhat_calib(eye_fid)

%Function to read only the calibration information from .eye files in the
%DOVES database

% Author    : Umesh Rajashekar, Ian van der Linde
% Date      : 27.07.01, revised 15.05.04


N = 9;

while(1)
    temp = fgets(eye_fid);
    if strfind(temp,'CALIBRATION POINTS')
        break;
    end
   
end

orig_calib_locs = zeros(N,2);
temp_calib = fscanf(eye_fid,'%f %f \n',N*2);	%read the x y and pupil of eye data
orig_calib_x = temp_calib(1:2:length(temp_calib));
orig_calib_y = temp_calib(2:2:length(temp_calib));
orig_calib_locs(:,1) = orig_calib_x;
orig_calib_locs(:,2) = orig_calib_y;

record_calib_locs = [];

while(~feof(eye_fid))
    
    temp = fgets(eye_fid);
    if strfind(temp,'Recalibration Data')
        temp_calib = fscanf(eye_fid,'%f %f \n',N*2);	%read the x y and pupil of eye data
        rec_calib_x = temp_calib(1:2:length(temp_calib));
        rec_calib_y = temp_calib(2:2:length(temp_calib));
        record_calib_locs = cat(3,record_calib_locs,[rec_calib_x rec_calib_y])   ;
    end
    
end



