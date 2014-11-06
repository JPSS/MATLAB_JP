%% startup
clc, clear all, close all
path0 = cd;
run('my_prefs')

%% choose colors
channel = cell(2,1);
channel{1} = 'red';
channel{2} = 'green';

%% LOAD STACK OF MOVIES
pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_ch1 = pickFirstFitsFiles(pname, channel{1}); 
files_ch2 = pickFirstFitsFiles(pname, channel{2});

N_movie = length(files_ch1);
if size(files_ch1,1) ~= size(files_ch2,1)
    disp('WARNING: not same number of movie files!')
end

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)


%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', ['Sequence ' channel{1} ':'], ['Sequence ' channel{2} ':'],... % sample options
    'Radius of peak [pixel]:', 'Number of frames for average [frames]:'};
input_default = {'2', '-1', '01', '10', '4', '10'};
tmp = inputdlg(input, 'Parameters', 1, input_default);
first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file
%determine sequences 
sequence_ch1 = zeros(1, size(tmp{3},2));
for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence_ch1(1,i) =1;
    end
end
sequence_ch2 = zeros(1, size(tmp{4},2));
for i=1:size(tmp{4},2)
    if(tmp{4}(i) == '1')
        sequence_ch2(1,i) =1;
    end
end
r_find = str2double(tmp(5)); % radius used to find spots
N_frames_for_average = str2double(tmp(6)); % minimal number of found spots in a trace

%% generate movie classes
ch1 = cell(N_movie,1);
ch2 = cell(N_movie,1);

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    ch2{i} = movie(pname, files_ch2{i}, first, last, sequence_ch2); % pname, fname, first, last, sequence
end

%% compute averages images
avg_img = cell(N_movie, 2);

for i=1:N_movie
    avg_img{i, 1} = ch1{i}.average_image(N_frames_for_average);
    avg_img{i, 2} = ch2{i}.average_image(N_frames_for_average);
end


%%



%% pick reference pattern
box_size = 16;
n_ref = 1;
mirror = 0;
img_size = size(avg_img{1,1}, 1);
dalpha = 5;


%% generate template


%{

templates = zeros(box_size_template, box_size_template, n_ref);
w = (box_size_template-1)/2;

close all
fig_dim =1.0*[20 20];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
imagesc(avg_img{1,1}), axis image, colormap gray

% select particle
h = imrect(gca, [img_size/2 img_size/2 double(box_size_template) double(box_size_template)]);
setResizable(h,0) 
pos = int16(wait(h));

%refine reference
r = double(box_size/2);
template = avg_img{1,1}(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3));

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
template = avg_img{1,1}(area(1):area(2), area(3):area(4));
%}



%% Calculate X-Correlation and find maximum correlations
disp('Calculating x-correlation...')
alpha = 0:dalpha:359;
n_rot = length(alpha); % number of rotations
dx = int16((box_size_template-box_size-1)/2);

peaks = cell(N_movie, 2);

% generate library
lib = zeros(box_size, box_size, length(alpha));
for j=1:n_rot

    lib(:,:,j) = generate_template(box_size, 4.8, 1.3, alpha(j)*pi/180); %box_size, distance between fluorophores, width of peak

    
    
end


%%
%loop through images
img3 = zeros(img_size, img_size, N_movie, 2); % stores maximum of cor-coef of all rotations
img3_index = zeros(img_size, img_size, N_movie, 2); % stores index of maximum

for i=1:N_movie
    for j=1:2 % loop thoru channels
        % channel j
        xcor_img = zeros(img_size, img_size, n_rot);

        for r=1:n_rot % loop through rotations
            tmp = normxcorr2(lib(:,:,r), avg_img{i,j}); % x-correlate
            xcor_img(:,:,r) = tmp(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);
        end
        % generate img3
        for k=1:512
        for l=1:512
            [cmax, imax] = max(xcor_img(k,l,:));
            img3(k,l,i,j) = cmax;
            img3_index(k,l,i,j) = imax;
        end
        end

        %find peaks in img3
        tmp = img3(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2,i,j);
        h_min = 0.7;

        p = find_peaks2d(tmp, round(box_size/4), h_min, 1); % radius of window, minimal height,  1= absolute
        if length(p) > 0
            p(:,1:2) = p(:,1:2)+box_size/2+1;
            tmp = img3(:,:,i,j);
            tmp_index = img3_index(:,:,i,j);
            idx = sub2ind(size(tmp), p(:,2), p(:,1) );
            peaks{i,j} = [p(:,1:2) tmp(idx) alpha(tmp_index(idx))' tmp_index(idx)]; % x y coer_coef alpha alpha_index

            display(['Movie (' num2str(i) ' of ' num2str(N_movie) '), channel ' channel{j} ' found ' num2str(size(p,1)) ' spots' ])

        end
    end

end
close all


%% generate stack of particles
close all
particles = cell(N_movie, 2);
particles_rot = cell(N_movie, 2);

for i=1:N_movie
    for j=1:2
    N = img_size;
    w = box_size/2;
    
    p = zeros(2*w+1, 2*w+1, size(peaks{i,j},1));
    p_rot = zeros(2*w+1, 2*w+1, size(peaks{i,j},1));

        for n=1:size(peaks{i,j},1) % loop through peak of image i,j
            x = peaks{i,j}(n,1);
            y = peaks{i,j}(n,2);
            angle = peaks{i,j}(n,4);
            p(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = avg_img{i,j}(max(1, y-w):min(N,y+w) , max(1, x-w):min(N,x+w) );
            
                %    lib(:,:,j) = tmp(dx:dx+box_size, dx:dx+box_size);
            
            w2=w+dx;
            
           % bla = zeros(2*w2+1, 2*w2+1);
          %  bla( max(1,w2-y+2):min(2*w2+1, 2*w2+1-y-w2+N),max(1,w2-x+2):min(2*w2+1, 2*w2+1-x-w2+N)) = img(max(1, y-w2):min(N,y+w2) , max(1, x-w2):min(N,x+w2) );
            tmp = imrotate(avg_img{i,j}(max(1, y-w2):min(N,y+w2) , max(1, x-w2):min(N,x+w2) ), -angle, 'crop');
          %  tmp = imrotate(bla, -angle, 'crop');

          %  p_rot(:,:,j) = bla(4*dx:4*dx+4*box_size, 4*dx:4*dx+4*box_size);
       
            p_rot(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = tmp(dx:dx+box_size, dx:dx+box_size);
            
            %{
            subplot(3,2,1:4)
            imagesc(avg_img{i,j}), colorbar, colormap gray, axis image, hold on
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
            

            subplot(1,2,1)
            imagesc(p(:,:,j)), colorbar, colormap gray, axis image
            title(['i=' num2str(i) ' j=' num2str(j) ' cc=' num2str(peaks{i,j}(n,3))   ])
            subplot(1,2,2)
            imagesc(p_rot(:,:,j)), colorbar, colormap gray, axis image
            title(['Particle ' num2str(j) ' angle = ' num2str(angle)])
            pause
            
            
            
        end
        particles{i,j} = p;
        particles_rot{i,j} = p_rot;

    end
end
clear('p', 'p_rot')



%%








