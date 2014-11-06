clear all, close all, clc
run('my_prefs'); path0=cd;
addpath([matlab_dir filesep 'TOOLS'])
addpath([matlab_dir filesep 'TEM_TOOLS' filesep 'spider_matlab'])

%% parameters
box_size = 50; % size of particle, HAS TO BE EVEN
dalpha = 5; %deg, angle resolution
mirror = 0; %include mirror transformation
img_size = 512;
%% load images
pname=uigetdir(data_dir,'Choose a folder with tem images.'); % get pathname
tmp = dir([pname filesep '*_16.TIF']);
fnames = {tmp.name}; % list of filenames
n_img = size(fnames,2);
path_out = [pname filesep 'particles']; % output folder
mkdir(path_out)

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
    h = imrect(gca, [size(img,1)/2 size(img,1)/2 double(box_size_template) double(box_size_template)]);
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
%%
close all
dx = (box_size_template-box_size-1)/2;
for i=1:n_ref
    subplot(n_ref,mirror+1,i)
    imagesc(templates(dx:dx+box_size, dx:dx+box_size,i)),  colormap gray, axis image
end



%% load images
images = zeros(512, 512, n_img);
for i=1:n_img
    img = imread([pname filesep fnames{i}]);
    images(:,:,i) = imresize(img(1:2048,1:2048),[512 512], 'nearest'); %bin image 4x4 for faster image processing
end

%% 
alpha = 0:dalpha:359;
n_rot = length(alpha); % number of rotations
dx = (box_size_template-box_size-1)/2;

peaks = cell(n_ref, n_img);
peaks2 = cell(n_ref, n_img);

h = waitbar(0,'Calculating x-correlation... ? time remaining');

for t=1:n_ref
    % generate library
    lib = zeros(box_size+1, box_size+1, length(alpha));
    for j=1:n_rot
        tmp = imrotate(templates(:,:,t), alpha(j), 'crop');
        lib(:,:,j) = tmp(dx:dx+box_size, dx:dx+box_size);
    end
 
    %loop through images
    img2 = zeros(img_size+box_size, img_size+box_size, n_img); % stores sum of quadratic xcoef
   % img2 = zeros(img_size, img_size, n_img); % stores sum of quadratic xcoef
    img3 = zeros(img_size, img_size, n_img); % stores maximum of cor-coef of all rotations
    img3_index = zeros(img_size, img_size, n_img); % stores index of maximum
    
    for i=1:n_img
        tic
      %  disp(['Reference ' num2str(t) ', image ' num2str(i)])
        xcor_img = zeros(512, 512, n_rot);
        
        for r=1:n_rot % loop through rotations
            tmp = normxcorr2(lib(:,:,r), images(:,:,i)); % x-correlate
            xcor_img(:,:,r) = tmp(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);
           % img2(:,:,i) = img2(:,:,i) + tmp.^2;
        end
        
        
        for k=1:512
        for l=1:512
            [cmax, imax] = max(xcor_img(k,l,:));
            img3(k,l,i) = cmax;
            img3_index(k,l,i) = imax;
        end
        end
        
        
        %find peaks in img3
        tmp = img3(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2,i);
        h_min = mean(tmp(:)) + 0.25*std(tmp(:));
        p = find_peaks(tmp, round(box_size/4), h_min, 0); % radius of window, minimal height,  no refinement
        p(:,1:2) = p(:,1:2)+box_size/2+1;
        
        tmp = img3(:,:,i);
        tmp_index = img3_index(:,:,i);
        idx = sub2ind(size(tmp), p(:,2), p(:,1) );
        peaks{t,i} = [p(:,1:2) tmp(idx) alpha(tmp_index(idx))']; % x y coer_coef alpha 
    
        display(['Reference ' num2str(t) ', image ' num2str(i) ', found ' num2str(size(p,1)) ' particles' ])

        if t==1 && i==1
            dt = toc;
        else
            dt_this = toc;
            dt = (dt+dt_this)/2;
        end
        frac = ((t-1)*n_img+i) / (n_ref*n_img);
        n_remain = (n_ref*n_img)-((t-1)*n_img+i);
        waitbar( frac , h, ['Calculating x-correlation... ' num2str(round(n_remain*dt/60*10)/10) 'min remaining'])
        
    end
    
    
end
close(h)

        
%%

%{

r = 1;
for i=1:n_img
tmp = img2(:,:,i);
h_min = mean(tmp(:)) + 0.25*std(tmp(:));
p1 = find_peaks(img2(:,:,i), round(box_size/4), h_min, 0); % radius of window, minimal height,  no refinement
tmp = img3(:,:,i);
h_min = mean(tmp(:)) + 0.25*std(tmp(:));
p2 = find_peaks(img3(:,:,i), round(box_size/4), h_min, 0); % radius of window, minimal height,  no refinement


close all
subplot(1, 2, 1)
%imagesc(img2(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2
imagesc(images(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2

plot(p1(:,1)+1-box_size/2,p1(:,2)+1-box_size/2,  'r.')

subplot(1, 2, 2)
%imagesc(img3(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2
imagesc(images(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2

plot(p2(:,1)+1,p2(:,2)+1,  'r.')
%}

%% view images and found particles

cc = varycolor(n_ref);
close all

for i=1:n_img
    subplot(1, 2, 1)
    imagesc(img3(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2
    for r=1:n_ref
        plot(peaks{r,i}(:,1), peaks{r,i}(:,2),  '.', 'color', cc(r,:))
    end
    title(num2str(i))
    hold off
    
    subplot(1, 2, 2)
    imagesc(images(:,:,i)), colorbar, colormap gray, axis image, hold on %
    for r=1:n_ref
        plot(peaks{r,i}(:,1), peaks{r,i}(:,2),  '.', 'color', cc(r,:))
    end
    title(num2str(i))
    hold off
    pause
end

            
%% Refine found particles
        

close all

for i=1:n_img
    tmp = imread([pname filesep fnames{i}]);
    img = imresize(tmp(1:2048,1:2048),[512 512], 'nearest'); %bin image 4x4 for faster image processing    
    N = 512;
    w = box_size/2;
    p= zeros(2*w+1, 2*w+1, size(peaks{t,i},1));
    for j=1:size(peaks{t,i},1)
        x = peaks{t,i}(j,1);
        y = peaks{t,i}(j,2);
        angle = peaks{t,i}(j,4);
        p(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = imrotate(img(max(1, y-w):min(N,y+w) , max(1, x-w):min(N,x+w) ), -angle, 'crop');
        
        subplot(1,2,1)
        imagesc(img), colorbar, colormap gray, axis image, hold on
        plot(x,y, 'r.')
        hold off
        subplot(1,2,2)
        imagesc(p(:,:,j)), colorbar, colormap gray, axis image
        title(['Particle ' num2str(j) ' angle = ' num2str(angle)])
        pause
        
    end

end


%%
t = 1; 
i = 1;
    tmp = imread([pname filesep fnames{i}]);
    img = imresize(tmp(1:2048,1:2048),[512 512], 'nearest'); %bin image 4x4 for faster image processing    
    N = 512;
    w = box_size/2;
    p= zeros(2*w+1, 2*w+1, 1, size(peaks{t,i},1));
for j=1:size(peaks{t,i},1)
        x = peaks{t,i}(j,1);
        y = peaks{t,i}(j,2);
        angle = peaks{t,i}(j,4);
        p(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),1,j) = imrotate(img(max(1, y-w):min(N,y+w) , max(1, x-w):min(N,x+w) ), -angle, 'crop');
        
end

close all
h = montage(uint16(p));
bla = ceil(sqrt(size(p,4)));
for  j=0:size(p,4)-1
    dx = 50;
    dy = 50;
    text(dx*(mod(j,bla))+25, floor(j/bla)*dy+45 , [ num2str(round(100*peaks{t,i}(j+1,3))/100) ], 'Color', [1 0 0])
end

%%
close all
particles = cell(n_ref, n_img);

for i=1:n_img
    tmp = imread([pname filesep fnames{i}]);
    %img = imresize(tmp(1:2048,1:2048),[512 512], 'nearest'); %bin image 4x4 for faster image processing    
    img = tmp(1:2048,1:2048);
    N = 2048;
    w = 4*box_size/2;
    %box_size*4
    for t=1:n_ref
        p= zeros(2*w+1, 2*w+1, size(peaks{t,i},1));
        for j=1:size(peaks{t,i},1)
            x = 4*peaks{t,i}(j,1);
            y = 4*peaks{t,i}(j,2);
            angle = peaks{t,i}(j,4);
            p(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = imrotate(img(max(1, y-w):min(N,y+w) , max(1, x-w):min(N,x+w) ), -angle, 'crop');
            %{
            subplot(1,2,1)
            imagesc(img), colorbar, colormap gray, axis image, hold on
            plot(x,y, 'r.')
            hold off
            subplot(1,2,2)
            imagesc(p(:,:,j)), colorbar, colormap gray, axis image
            title(['Particle ' num2str(j) ' angle = ' num2str(angle)])
            pause
            %}
        end
        particles{t,i} = p;
    end
end
%% display images
%{
close all
for t=1:n_ref
    for i=1:n_img
        for j=1:size(particles{t,i},3)
            imagesc(particles{t,i}(:,:,j)), colorbar, colormap gray, axis image
            title(['Image ' num2str(i) ' particle ' num2str(j)])
            pause
        end
    end
end
%} 

%% write particles for each reference
for t=1:n_ref
    n_particle = 0;
    for i=1:n_img
        n_particle = n_particle + size(particles{t,i}, 3);
    end    
    p_out = zeros(2*w+1, 2*w+1, n_particle);

    m=1;
    for i=1:n_img
        for j=1:size(particles{t,i},3)
            p_out(:,:,m) = particles{t,i}(:,:,j);
            m = m+1;
        end
    end
    writeSPIDERfile([path_out filesep 'ref_' num2str(t) '.spi'], p_out, 'stack')

end

%%
disp('finished')
