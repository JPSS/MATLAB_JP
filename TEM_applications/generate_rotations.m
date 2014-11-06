%% statrup
close all, clear all, clc
run('my_prefs')
path0 = cd;

%% get file
cd(data_dir)
[fname pname] = uigetfile('.png', 'Select particle to be randomly rotated');
cd(path0)

img = imread([pname fname]);

path_out = [pname fname(1:end-4) '_rotated' filesep];
mkdir(path_out);

%% generate roation of this particle
N = 100;

angle = 360*rand(N,1); % random angles between 0 and 360
flip = round(rand(N,1));

stack = uint32(zeros(size(img,1), size(img,2), N));
for i=1:N
    if flip(i) == 1
        tmp = flipdim(img,1);
    else
        tmp = img;
    end
    stack(:,:,i) = imrotate(tmp, angle(i), 'crop');
end

%% write output
for i=1:N
    imwrite(uint16(stack(:,:,i)), [path_out fname(1:end-4) '_' sprintf('%.2i.png', i)])
end
disp('Done')