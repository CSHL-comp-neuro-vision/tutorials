%Example to read calibration data directly from the .eye file in DOVES

% Author    : Umesh Rajashekar, Ian van der Linde
% Date      : 27.07.01, revised 15.05.04s


clear;close all;

% Define Constants
IMAGES_DIR      = '../Images/';
DATA_DIR        = '../RawData/';


%% Select a subject's .eye file
[file_name DATA_DIR] = uigetfile([DATA_DIR '*.eye'],'Please select a Subject');


%% Process Calibration information
eye_fid = fopen([DATA_DIR file_name],'r');
[orig_calib_locs, record_calib_locs]  = read_vanhat_calib(eye_fid);

figure;
%Plot the actual location of the calibration dots
plot(orig_calib_locs(:,1), orig_calib_locs(:,2), 'ro', 'MarkerSize', 16,'MarkerFaceColor','r');
set(gca,'FontSize',14);
hold on; axis ij; axis([1 1024 1 768]);
title(['Calibration for ' file_name]);
for iRecal = 1:size(record_calib_locs,3)
    %Plot the corresponding location of the observer's fixations
    h = plot(record_calib_locs(:,1, iRecal), record_calib_locs(:,2, iRecal), 'ko', 'MarkerSize', 8,'MarkerFaceColor','k');
    legend('Actual location', 'Recorded location');
    display('Hit any key to continue');
    pause;
    delete(h);
end

fclose(eye_fid);

%% 