%% startup
clear all; close all; clc;
run('my_prefs')

%% select folder
pname = uigetdir(data_dir, 'Select directory with .xmp files');
files = dir([pname filesep 'ml2d_ref*.xmp']); % get files in folder
path_out = [pname filesep 'averages_tif']; % generate output folder
mkdir(path_out)

%% save as tif images
for i=1:size(files,1)
    img = readSPIDERfile([pname filesep files(i).name]);
    img = img - min(img(:));
    imwrite(uint16(img*(2^16-1)./max(img(:))), [path_out filesep files(i).name(1:end-4) '.tif']);
end
    
disp(['Done... converted ' num2str(size(files,1)) ' files.'])
