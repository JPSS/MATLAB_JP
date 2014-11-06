%% startup
close all, clear all, clc
run('my_prefs')

%% get pathname 
pname = uigetdir(data_dir, 'Select folder, where green/red/blue is stored');

colors = ones(1,3); % r g b

dir_red = dir([pname filesep 'red*']);
if isempty(dir_red)
    colors(1) = 0;
else
    path_red = [pname filesep dir_red.name];
    disp(['Found red folder: ' path_red])
end

dir_green = dir([pname filesep 'green*']);
if isempty(dir_green)
    colors(2) = 0;
else
    path_green = [pname filesep dir_green.name];
    disp(['Found green folder: ' path_green])
end

dir_blue = dir([pname filesep 'blue*']);
if isempty(dir_blue)
    colors(3) = 0;
else
    path_blue = [pname filesep dir_blue.name];
    disp(['Found blue folder: ' path_blue])
end

%% get stem folders from red channel

tmp = dir(path_red);
dir_red = cell(0,1); % list of folders in path_red
for i=1:length(tmp)
    if tmp(i).isdir && ~strcmp(tmp(i).name, '.') && ~strcmp(tmp(i).name, '..')
        dir_red = [dir_red; tmp(i).name];
    end
end

%%  move files
for i=1:length(dir_red)
    path_out = [pname filesep dir_red{i}];
    mkdir(path_out); 
    
    % move red files
    if colors(1)
        path_source = [path_red filesep dir_red{i} ];
        if exist(path_source, 'dir') == 7 % directory exists
            movefile([path_source filesep '*'], path_out)
            rmdir(path_source); % remove source directory
        end
        
    end
    
    % move green files
    if colors(2)
        path_source = [path_green filesep dir_red{i} ];
        if exist(path_source, 'dir') == 7 % directory exists
            movefile([path_source filesep '*'], path_out)
            rmdir(path_source); % remove source directory
        end
        
    end
    
    % move blue files
    if colors(3)
        path_source = [path_blue filesep dir_red{i} ];
        if exist(path_source, 'dir') == 7 % directory exists
            movefile([path_source filesep '*'], path_out)
            rmdir(path_source); % remove source directory
        end       
    end
    
end

disp('Done moving files')
        
