clear all, close all, clc
run('my_prefs'); path0=cd;
%% parameters
box_size = 50; % size of particle, HAS TO BE EVEN
dalpha = 5; %deg, angle resolution
mirror = 0; %include mirror transformation

%% load images
pname=uigetdir(data_dir,'Choose a folder with tem images.'); % get pathname
tmp = dir([pname filesep '*_16.TIF']);
fnames = {tmp.name}; % list of filenames
n_img = size(fnames,2);
%% create class references
tmp = inputdlg({'Number of class references:'}, 'Number of classes', 1, {'1'});
n_ref = str2double(tmp(1));

box_size_template = ceil(2*box_size/2/cos(pi/4));%200;
box_size_template = int16(box_size_template + mod(box_size_template, 2) +1 ) ; % make it even
templates = zeros(box_size_template, box_size_template, n_ref);
w = (box_size_template-1)/2;
for i=1:n_ref 
    go_on = 1;
    j = 1;
    close all
    fig_dim =1.0*[10 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

    while go_on 
        img = imread([pname filesep fnames{j}]);
        img = img(1:2048,1:2048);
        img = imresize(img,[512 512], 'nearest'); %bin image 4x4 for faster image processing
        imagesc(img), axis image, colormap gray
        button = questdlg(['Select an image for ref. ' num2str(i)],'Image','Use this','Previous','Next', 'Use this');
        if strcmp(button, 'Next')
            j = min(n_img,j+1);
        end
        if strcmp(button, 'Previous')
            j = max(j-1, 1);
        end
        if strcmp(button, 'Use this')
            go_on = 0;
        end
    end
    % select particle
    h = imrect(gca, [img_size(1)/2 img_size(1)/2 double(box_size_template) double(box_size_template)]);
    setResizable(h,0) 
    pos = int16(wait(h));

    % refine reference
    close all
    template = img(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3));
    imagesc( [pos(1) pos(1)+pos(3)], [pos(2) pos(2)+pos(4)],template), colorbar, colormap gray, axis image
    c = round(ginput(1));
    area = [c(2)-w c(2)+w c(1)-w c(1)+w];
    templates(:,:,i) = img(area(1):area(2), area(3):area(4));

  %  imagesc([area(1) area(2)], [area(3) area(4)],template), colorbar, colormap gray, axis image, hold on
  %  plot(c(2), c(1), 'r.')  
end

%% generate template library
alpha = 0:dalpha:359;
N_templates = length(alpha);
dx = (box_size_template-box_size-1)/2;
lib = zeros(box_size+1, box_size+1, length(alpha), n_ref);
for i = 1:n_ref
    
    for j=1:length(alpha)
         tmp = imrotate(templates(:,:,i), alpha(j), 'crop');
         lib(:,:,j,i) = tmp(dx:dx+box_size, dx:dx+box_size);
    end
end

%% view some templates
%{
close all
for i=1:n_ref
    for j=1:length(alpha)
            imagesc(lib(:,:,j,1)), axis image, colormap gray
            title(['Reference ' num2str(i) ', alpha = ' num2str(alpha(j))])
            pause
    end
end
%}
%% calculate corss-correlation matrices
ximg = zeros(img_size(1)+box_size, img_size(2)+box_size, N_templates, n_ref);

for j=1:n_ref
    h = waitbar(0,['Calc. xcorr for ref ' num2str(j)]); s = clock;
    for i = 1:N_templates
        ximg(:,:,i,j) = normxcorr2(lib(:,:,i,j), img); % calculate xcorrelation function

        if i ==1 % begin estimate remaining time
            is = etime(clock,s);
            esttime = is * N_templates;
        end

        if i == 10 % begin estimate remaining time
            is = etime(clock,s)/10;
            esttime = is * N_templates;
        end
        h = waitbar(i/N_templates,h,['Correlating ref. ' num2str(j) ', time remaining = ' num2str(esttime-etime(clock,s),'%4.1f') ' sec' ]); %update waitbar
    end
    close(h)
end

%% create average and squared image
img2 = zeros(size(ximg,1), size(ximg,2), n_ref);
for i=1:n_ref
    for j=1:N_templates
       img2(:,:,i) = img2(:,:,i) + ximg(:,:,j, i).^2; 
    end
end

%% find peaks in img2
h_min = mean(img2(:)) + 0.5*std(img2(:));
display('Searching for peaks ...')
p = find_peaks(img2(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2), round(box_size/2), h_min, 0); % radius of window, minimal height,  no refinement
display(['found ' num2str(size(p,1)) ' peaks.'])

p_xy = p(:,1:2)+1%;-box_size/2;

%% plot result
close all
subplot(1, 2, 1)
imagesc(sqrt(img2)), colorbar,  colormap gray, axis image, hold on % img2
for j=1:size(p,1)
   plot(p(j,1)+1, p(j,2)+1, 'ro')
end

subplot(1, 2, 2)
imagesc(img), colorbar, colormap gray, axis image, hold on %
for j=1:size(p,1)
    plot(p_xy(j,1), p_xy(j,2), 'r.')
end

%% generate stack of particles
N_particles = size(p_xy,1);
picked = zeros(box_size+1, box_size+1, N_particles);

close all
for i=1:N_particles
    if p_xy(i,1)>0 && p_xy(i,1)<=512 && p_xy(i,2)>0 && p_xy(i,2)<=512
        index = [max(p_xy(i,2)-box_size/2, 1) min(p_xy(i,2)+box_size/2, 512) max(p_xy(i,1)-box_size/2, 1) min(p_xy(i,1)+box_size/2, 512)] ;
        tmp = img(index(1):index(2), index(3):index(4));
        a=ones(1,4);
        height = index(2)-index(1);
        
        width = index(4)-index(3);

        if width < 50
            if index(1) == 1 %paritlce at left corner
                a(1) = box_size-index(2)+1;
                a(2) = box_size;
            else  %paritlce at right corner
                a(1) = 1;
                a(2) = width+1;
            end
        else
            a(1) = 1;
            a(2) = box_size+1;
        end
        if height < 50
            if index(3) == 1 %paritlce at left corner
                a(3) = box_size-index(4)+1;
                a(4) = box_size;
            else  %paritlce at right corner
                a(3) = 1;
                a(4) = height+1;
            end
        else
            a(3) = 1;
            a(4) = box_size+1;
        end
            
        
        picked( a(3):a(4), a(1):a(2),i) = tmp;
            
        
        imagesc(picked(:,:,i)), colormap gray, axis image
        pause
    end
end


%% writing output



writeSPIDERfile('/Users/jonasfunke/Documents/MATLAB/TEM tools/mystack.spi', picked, 'stack')


