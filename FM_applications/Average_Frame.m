%% STARTUP
clc, clear all, close all
path0 = cd;
run('my_prefs')


%% select file
cd(data_dir)
[fname pname]=uigetfile('*.fits','Select a movie (.fits).');
cd(path0)

path_out = pname ;
prefix = fname(1:end-5);
%% generate a movie file
input = {'First Frame:', 'Last Frame (-1=all):', 'Sequence:', 'First n frames for average [frames]:'};
input_default = {'2', '-1', '10', '-1'};
tmp = inputdlg(input, 'Parameters', 1, input_default);

first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file

%determine sequences 
sequence = zeros(1, size(tmp{3},2));

for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence(1,i) =1;
    end
end

N_frames = str2double(tmp(4)); % minimal number of found spots in a trace
mov = movie(pname, fname, first, last, sequence); % pname, fname, first=1, last=all, sequence=all


%% generate average frame
avg_img = mov.average_image(N_frames);
avg_img_all = mov.average_image(-1);

%% write averages images
 
imwrite(uint16(avg_img), [path_out filesep prefix '_avg_' num2str(N_frames) '-frames.tif']);
imwrite(uint16(avg_img_all), [path_out filesep prefix '_avg_all-frames.tif']);


disp('Done.')