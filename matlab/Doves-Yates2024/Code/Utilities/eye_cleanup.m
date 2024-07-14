function[x_data,y_data,pupil_data] =  eye_cleanup(x_data,y_data,pupil_data,SCREEN_WIDTH,SCREEN_HEIGHT)

%First clean up the data
x_index = find( (x_data > SCREEN_WIDTH) | (x_data < 1));
x_data(x_index)=[];y_data(x_index) = [];pupil_data(x_index) = [];

y_index = find((y_data>SCREEN_HEIGHT) | (y_data <= 0));
x_data(y_index)=[];y_data(y_index) = [];pupil_data(y_index) = [];

pupil_index = find(pupil_data == 0);
x_data(pupil_index) =[];y_data(pupil_index) =[];pupil_data(pupil_index) =[];
