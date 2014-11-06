%% startup
clear all; close all; clc;
run('my_prefs'); path0=cd;

%% parameters
input = {'Box size [pixel]:', 'Number of class references:', 'Include mirror (0=no, 1=yes):',... % sample options
    'Binning:', 'Radius of high-pass filter [pixel]:', 'Angel resolution [deg]:'};
input_default = {'200', '1', '1', '4', '15', '5'};
tmp = inputdlg(input, 'Parameters', 1, input_default);
box_size_real = round(str2double(tmp(1))); % size of particle on real image, HAS TO BE EVEN
mirror = str2double(tmp(3)); % include mirror transformation, 0=no, 1= yes
n_ref = str2double(tmp(2))*(mirror+1);
n_bin = str2double(tmp(4)); % number of pixel to bin in one dim
r_filter = str2double(tmp(5));  %100; % pixel (original image), radius for gaussian high pass, -1 == no filtering
dalpha = str2double(tmp(6)); % deg, angle resolution

if mod(box_size_real,n_bin*2) ~= 0
    box_size_real = box_size_real+2*n_bin-mod(box_size_real,n_bin*2);
    disp(['Box size not a multiple of bin_size. Changed to ' num2str(box_size_real)])
end

box_size = box_size_real/n_bin; % size of particle, HAS TO BE EVEN
img_size = 2048/n_bin; %size of the binned image


%% load images
pname=uigetdir(data_dir,'Choose a folder with tem images.'); % get pathname
tmp = dir([pname filesep '*_16.TIF']);
fnames = {tmp.name}; % list of filenames
n_img = size(fnames,2);
path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_particles']; % output folder
mkdir(path_out)


%% load images
disp(['Loading and filtering ' num2str(n_img) ' images...'])
images = zeros(img_size, img_size, n_img);
if r_filter > 0
    f_filter = fspecial('gaussian', r_filter*4*2 , r_filter); % gaussian filter, diameter = 2*(width = 4*sigma)
end
h = waitbar(0,'Loading and filtering images... ? time remaining');
tic
for i=1:n_img
    img = imread([pname filesep fnames{i}], 'PixelRegion', {[1 2048], [1 2048]});
    if r_filter > 0
        tmp = double(img)-double(imfilter(img, f_filter, 'same'));
    else
        tmp = double(img);
    end
    images(:,:,i) = imresize(tmp,[img_size img_size], 'nearest'); %bin image 4x4 for faster image processing
    %images(:,:,i) = imresize(img(1:2048,1:2048),[512 512], 'nearest'); %bin image 4x4 for faster image processing
    if i==1
        dt = toc;
        tic;
    else
        dt_this = toc;
        dt = (dt+dt_this)/2;
        tic;
    end
    
    n_remain = n_img-i;
    waitbar( i/n_img , h, ['Loading and filtering images... ' num2str(round(n_remain*dt/60*10)/10) ' min remaining'])

end
toc;
close(h); close all;
%% create class references

box_size_template = ceil(2*box_size/2/cos(pi/4));%200;
box_size_template = int16(box_size_template + mod(box_size_template, 2) +1 ) ; % make it even
templates = zeros(box_size_template, box_size_template, n_ref);
w = (box_size_template-1)/2;
for i=1:n_ref/(mirror+1) 
    go_on = 1;
    j = 1;
    close all
    fig_dim =1.0*[10 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

    while go_on 
        imagesc(images(:,:,j)), axis image, colormap gray
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
    h = imrect(gca, [img_size/2 img_size/2 double(box_size_template) double(box_size_template)]);
    setResizable(h,0) 
    pos = int16(wait(h));

    
    %refine reference
    r = double(box_size/2);
    template = images(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),j);

    close all
    imagesc( [pos(1) pos(1)+pos(3)], [pos(2) pos(2)+pos(4)],template), colorbar, colormap gray, axis image
    cur_fig = gca;
    h = impoint(gca,[pos(1)+box_size_template/2 pos(2)+box_size_template/2]);
    setColor(h, [1 0 0])
    b = ellipse(r, r, 0, double(pos(1)+box_size_template/2), double(pos(2)+box_size_template/2));
    set(b, 'Color', [1 0 0])
    addNewPositionCallback(h, @(pos) update_ellipse(pos, r, cur_fig) );
    c = round(wait(h));
    close all
    

    area = [c(2)-w c(2)+w c(1)-w c(1)+w];
    if mirror
        templates(:,:,2*i-1) = images(area(1):area(2), area(3):area(4), j);
        templates(:,:,2*i) = flipdim(images(area(1):area(2), area(3):area(4), j) ,1);
    else
        templates(:,:,i) = images(area(1):area(2), area(3):area(4), j);
    end
      
end
close all

%% display and write templates
path_out_templates = [path_out filesep 'reference_particles'];
mkdir(path_out_templates)
dx = (box_size_template-box_size-1)/2;
close all
for i=1:n_ref
    subplot(n_ref/(mirror+1),mirror+1,i)
    imagesc(templates(dx:dx+box_size, dx:dx+box_size,i)),  colormap gray, axis image
    tmp_out = templates(dx:dx+box_size, dx:dx+box_size,i)-min(min(templates(dx:dx+box_size, dx:dx+box_size,i)));
    imwrite(  uint16(tmp_out*(2^16-1)/max(tmp_out(:))) , [path_out_templates filesep 'ref_' num2str(i) '.tif' ]);
end

%% ask if one wants to refine selection
refine = 1; %strcmp(questdlg('Dismiss particles with low correlation?','Refinement','Yes','No','No'), 'Yes');  

%% Calculate X-Correlation and find maximum correlations
disp('Calculating x-correlation...')
alpha = 0:dalpha:359;
n_rot = length(alpha); % number of rotations
dx = int16((box_size_template-box_size-1)/2);

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
    img3 = zeros(img_size, img_size, n_img); % stores maximum of cor-coef of all rotations
    img3_index = zeros(img_size, img_size, n_img); % stores index of maximum
    
    for i=1:n_img
        tic
        xcor_img = zeros(img_size, img_size, n_rot);
        
        for r=1:n_rot % loop through rotations
            tmp = normxcorr2(lib(:,:,r), images(:,:,i)); % x-correlate
            xcor_img(:,:,r) = tmp(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);
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
    
        p = find_peaks2d(tmp, round(box_size/4), h_min, 0); % radius of window, minimal height,  no absolure = relative height
        if length(p) > 0
            p(:,1:2) = p(:,1:2)+box_size/2+1;
            tmp = img3(:,:,i);
            tmp_index = img3_index(:,:,i);
            idx = sub2ind(size(tmp), p(:,2), p(:,1) );
            peaks{t,i} = [p(:,1:2) tmp(idx) alpha(tmp_index(idx))']; % x y coer_coef alpha 

            display(['Reference (' num2str(t) '/' num2str(n_ref) '): image (' num2str(i) '/' num2str(n_img) '), found ' num2str(size(p,1)) ' particles' ])

        end

        if t==1 && i==1
            dt = toc;
        else
            dt_this = toc;
            dt = (dt+dt_this)/2;
        end
        frac = ((t-1)*n_img+i) / (n_ref*n_img);
        n_remain = (n_ref*n_img)-((t-1)*n_img+i);
        waitbar( frac , h, ['Calculating x-correlation... ' num2str(round(n_remain*dt/60*10)/10) ' min remaining']) 
    end
end
pause(0.1)
close(h)
close all














%% remove particles, which belong to multiple classes
%cc = varycolor(n_ref);
peaks_ref = cell(n_ref, n_img);
h = waitbar(0,'Searching for particles... ');

for i=1:n_img
    disp(['refining ' num2str(i) ])
    
    % genertate image of  correlations
    cor_img = zeros(img_size,img_size );
    cor_img_index = zeros(img_size,img_size );
    rot_img = zeros(img_size,img_size );
    for r=1:n_ref
        idx = sub2ind(size(cor_img), peaks{r,i}(:,2), peaks{r,i}(:,1) );
        cor_img(idx) = peaks{r,i}(:,3);
        cor_img_index(idx) = r;
        rot_img(idx) = peaks{r,i}(:,4);

    end
    
    
    
    p = find_peaks2d(cor_img, round(box_size/4), 0, 1 ); % find-peaks, width, min_height, absolute height 
    p(:,1:2) =  p(:,1:2)+1;
   
    
    for j=1:size(p,1)
        peaks_ref{cor_img_index(p(j,2),p(j,1)),i} = [peaks_ref{cor_img_index(p(j,2),p(j,1)),i}; p(j,1:2) cor_img(p(j,2),p(j,1)) rot_img(p(j,2),p(j,1)) ];
    end
    
    %{
     close all
    imagesc(cor_img), colorbar,  colormap gray, axis image, hold on % img2

    for r=1:n_ref
        plot(peaks{r,i}(:,1), peaks{r,i}(:,2),  'o', 'color', cc(r,:))
        plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '.', 'color', cc(r,:));

    end
    pause 
    %}
    
    waitbar( i/n_img , h, 'Searching for particles... ')

    
end
close(h);
pause(0.1);
    
   

%% display images and found particles
%{
cc = varycolor(n_ref);
close all
myleg = cell(n_ref,1);
for t=1:n_ref
    myleg{t} = ['Reference ' num2str(t)];
end
for i=1:n_img
    subplot(1, 2, 1)
    imagesc(img3(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2
    for r=1:n_ref
        h(r) = plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '.', 'color', cc(r,:));
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
    
    subplot(1, 2, 2)
    imagesc(images(:,:,i)), colorbar, colormap gray, axis image, hold on %
    for r=1:n_ref
        h(r)=plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '.', 'color', cc(r,:));
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
    pause
end
%}

%% generate stack of particles
close all
particles = cell(n_ref, n_img);
particles_rot = cell(n_ref, n_img);

for i=1:n_img
    img = imread([pname filesep fnames{i}], 'PixelRegion', {[1 2048], [1 2048]});
    N = 2048;
    w = 4*box_size/2;
    for t=1:n_ref
        p= zeros(2*w+1, 2*w+1, size(peaks_ref{t,i},1), 'uint16');
        p_rot= zeros(2*w+1, 2*w+1, size(peaks_ref{t,i},1), 'uint16');

        for j=1:size(peaks_ref{t,i},1)
            x = 4*peaks_ref{t,i}(j,1);
            y = 4*peaks_ref{t,i}(j,2);
            angle = peaks_ref{t,i}(j,4);
            p(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = img(max(1, y-w):min(N,y+w) , max(1, x-w):min(N,x+w) );
            
                %    lib(:,:,j) = tmp(dx:dx+box_size, dx:dx+box_size);
            
            w2=w+4*dx;
            
           % bla = zeros(2*w2+1, 2*w2+1);
          %  bla( max(1,w2-y+2):min(2*w2+1, 2*w2+1-y-w2+N),max(1,w2-x+2):min(2*w2+1, 2*w2+1-x-w2+N)) = img(max(1, y-w2):min(N,y+w2) , max(1, x-w2):min(N,x+w2) );
            tmp = imrotate(img(max(1, y-w2):min(N,y+w2) , max(1, x-w2):min(N,x+w2) ), -angle, 'crop');
          %  tmp = imrotate(bla, -angle, 'crop');

          %  p_rot(:,:,j) = bla(4*dx:4*dx+4*box_size, 4*dx:4*dx+4*box_size);
       
            p_rot(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = tmp(4*dx:4*dx+4*box_size, 4*dx:4*dx+4*box_size);
            %{
            subplot(3,2,1:4)
            imagesc(img), colorbar, colormap gray, axis image, hold on
            plot(x,y, 'r.')
            hold off
            subplot(3,2,5)
            imagesc(p(:,:,j)), colorbar, colormap gray, axis image
            title(['Particle ' num2str(j) ' angle = 0' ])
            subplot(3,2,6)
            imagesc(p_rot(:,:,j)), colorbar, colormap gray, axis image
            title(['Particle ' num2str(j) ' angle = ' num2str(angle)])
            pause
            %}
            
        end
        particles{t,i} = p;
        particles_rot{t,i} = p_rot;

    end
end
clear('p', 'p_rot')

%%
clear data
data(n_ref/(mirror+1)) = struct('stats', [], 'particles', [],  'particles_rot', []);
w = 4*box_size/2;

for t=1:n_ref/(mirror+1)
    n_particle = 0;
    for i=1:n_img
        n_particle = n_particle + size(particles{(mirror+1)*t,i}, 3);
        if mirror
            n_particle = n_particle + size(particles{(mirror+1)*t-1,i}, 3);
        end
    end    
    
    p_out = zeros(2*w+1, 2*w+1, n_particle, 'uint16');
    p_out_rot = zeros(2*w+1, 2*w+1, n_particle, 'uint16');

    stats_out = zeros(n_particle, 5);
    
    m=1;
    if mirror 
        for i=1:n_img
            for j=1:size(particles{(mirror+1)*t-1,i},3)
                p_out(:,:,m) = particles{(mirror+1)*t-1,i}(:,:,j);
                p_out_rot(:,:,m) = particles_rot{(mirror+1)*t-1,i}(:,:,j);

                stats_out(m,1:4)  = peaks_ref{(mirror+1)*t-1,i}(j,1:4);
                stats_out(m,5)  = i;
                m = m+1;
            end
        end
    end
    for i=1:n_img
        for j=1:size(particles{(mirror+1)*t,i},3)
            p_out(:,:,m) = particles{(mirror+1)*t,i}(:,:,j);
            if mirror
                p_out_rot(:,:,m) = flipdim(particles_rot{(mirror+1)*t,i}(:,:,j),1);
            else
                p_out_rot(:,:,m) = particles_rot{(mirror+1)*t,i}(:,:,j);
            end
            stats_out(m,1:4)  = peaks_ref{(mirror+1)*t,i}(j,1:4);
            stats_out(m,5)  = i;
            m = m+1;
        end
    end
    
    data(t).stats = stats_out;
    data(t).particles = p_out;
    data(t).particles_rot = p_out_rot;

end
clear('p_out', 'p_out_rot')

%% refine 
disp('Refining particles...')

limit = zeros(n_ref/(mirror+1), 1);
if refine
    for t=1:n_ref/(mirror+1)
        [cc_sort, sort_index] = sortrows(data(t).stats(:,3), -1);
        i = [1:size(cc_sort,1)]';

        close all
        subplot(1, 2, 1)
        imagesc(data(t).particles_rot(:,:,sort_index(1))), axis image, colormap gray
       % title(['Ref ' num2str(t) ', particle ' num2str(i) ', cc = ' num2str(data(t).stats(sort_index(1),3))])
        cur_img = gca;

        subplot(1, 2, 2)
        plot(i, cc_sort, 'b'), hold on
        ylim = [0 1];
        set(gca, 'YLim', ylim, 'XLim', [1 i(end)]);
        h = imline(gca,[1 1], ylim);
        setColor(h,[1 0 0]);
        setPositionConstraintFcn(h, @(pos)[ min( i(end), max(1,[pos(2,1);pos(2,1)])) ylim'   ])

        id = addNewPositionCallback(h, @(pos) update_img(  data(t).particles_rot(:,:,sort_index(  max(1, min(i(end), round(pos(1,1))))   ) ), cur_img )  );
        id2 = addNewPositionCallback(h, @(pos) title(['cc = ' num2str( cc_sort(  max(1, min(i(end), round(pos(1,1))))) )]) );

        pos_line = wait(h);
        limit_index = max(1, min(i(end), round(pos_line(1,1))));
        limit(t) = cc_sort(limit_index);

    end
    pause(0.1)
    close all
end
pause(0.1)
close all

%% write  refining 
dlmwrite([path_out filesep 'cc_thresholds_ref_cc_min.txt'], [ [1:n_ref/(mirror+1)]' limit], '\t');

%%

data_refined(n_ref/(mirror+1)) = struct('stats', [], 'particles', [],  'particles_rot', []);

w = 4*box_size/2;
for t=1:n_ref/(mirror+1)
    
    index = find(data(t).stats(:,3)>limit(t));
    stats = zeros(size(index, 1),5);

    p = zeros(2*w+1, 2*w+1, size(index, 1), 'uint16');
    p_rot = zeros(2*w+1, 2*w+1, size(index, 1), 'uint16');

    for i=1:size(index,1)
        stats(i,:) = data(t).stats(index(i),:);
        p(:,:,i) = data(t).particles(:,:,index(i));
        p_rot(:,:,i) = data(t).particles_rot(:,:,index(i));
    end
    
    data_refined(t).stats = stats;
    data_refined(t).particles = p;
    data_refined(t).particles_rot = p_rot;
    
end
clear('p', 'p_rot')


%% write particles for each reference 
disp('Writing particles...')

for t=1:n_ref/(mirror+1)
    % write as spider-file
    writeSPIDERfile([path_out filesep 'ref_' num2str(t) '.spi'], data_refined(t).particles, 'stack')
    writeSPIDERfile([path_out filesep 'ref_' num2str(t) '_rot.spi'], data_refined(t).particles_rot, 'stack')
  
    % write as single tif-files
    path_out_tif1 = [path_out filesep 'ref_' num2str(t) '_tif'];
    path_out_tif2 = [path_out filesep 'ref_' num2str(t) '_rot_tif'];
    path_out_spi = [path_out filesep 'ref_' num2str(t) '_spi'];
    
    mkdir(path_out_tif1)
    mkdir(path_out_tif2)
    mkdir(path_out_spi)
    
    for i=1:size(data_refined(t).particles, 3)
        imwrite(data_refined(t).particles(:,:,i), [path_out_tif1 filesep 'ref_' num2str(t) '_' sprintf('%.3i',i) '.tif' ]);
        imwrite(data_refined(t).particles_rot(:,:,i), [path_out_tif2 filesep 'ref_' num2str(t) '_' sprintf('%.3i',i) '.tif' ]);
        writeSPIDERfile([path_out_spi filesep 'ref_' num2str(t) '_' sprintf('%.3i',i) '.spi' ], data_refined(t).particles(:,:,i))

    end
end

%% write all particles
% write as spider-file
writeSPIDERfile([path_out filesep 'all.spi'], cat(3,data.particles), 'stack')

% write as single tif-files
path_out_tif = [path_out filesep 'all_tif'];
mkdir(path_out_tif)
m = 1;
for t=1:n_ref/(mirror+1)
    for i=1:size(data(t).particles, 3)
        imwrite(data(t).particles(:,:,i), [path_out_tif filesep 'all_' sprintf('%.3i',m) '.tif' ]);
        m = m+1;
    end
end

%% display images and found particles
disp('Plotting images...')

cc = varycolor(n_ref);
close all
fig_dim =2*[20 10];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

myleg = cell(n_ref,1);
for t=1:n_ref
    myleg{t} = ['Ref. ' num2str(t)];
end

w = box_size/2;

for i=1:n_img
    subplot(1, 2, 1)
    imagesc(img3(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2
    h = zeros(n_ref,1);
    for r=1:n_ref
        h(r) = plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '+', 'color', cc(r,:), 'MarkerSize', 1);
        
        for j=1:size(peaks_ref{r,i},1)
            rectangle('Position',[peaks_ref{r,i}(j,1)-w, peaks_ref{r,i}(j,2)-w, 2*w+1, 2*w+1], 'EdgeColor', cc(r,:))
        end
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
    
    subplot(1, 2, 2)
    imagesc(images(:,:,i)), colorbar, colormap gray, axis image, hold on %
    for r=1:n_ref
        h(r)=plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  'x', 'color', cc(r,:), 'MarkerSize', 1);
        for j=1:size(peaks_ref{r,i},1)
            rectangle('Position',[peaks_ref{r,i}(j,1)-w, peaks_ref{r,i}(j,2)-w, 2*w+1, 2*w+1], 'EdgeColor', cc(r,:))
        end
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
       
    %pause
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image2_' sprintf('%.03i',i) '.tif'])

end


%%
disp('finished')

