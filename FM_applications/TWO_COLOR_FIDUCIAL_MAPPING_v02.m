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

%% generate library
box_size = 16;
n_ref = 1;
mirror = 0;
img_size = size(avg_img{1,1}, 1);
dalpha = 2;

alpha = 0:dalpha:359;
n_rot = length(alpha); % number of rotations

% generate library
lib_1 = zeros(box_size, box_size, length(alpha));
for j=1:n_rot
    lib_1(:,:,j) = generate_template(box_size, 4.8, 1.5, alpha(j)*pi/180); %box_size, distance between fluorophores, width of peak
end
lib_2 = zeros(box_size, box_size, length(alpha));
for j=1:n_rot
    lib_2(:,:,j) = generate_template(box_size, 4.8, 1.5, alpha(j)*pi/180+pi); %box_size, distance between fluorophores, width of peak
end

%% Calculate X-Correlation and find maximum correlations
%loop through images
img3 = zeros(img_size, img_size, N_movie, 2); % stores maximum of cor-coef of all rotations
img3_index = zeros(img_size, img_size, N_movie, 2); % stores index of maximum
peaks = cell(N_movie, 2);

for i=1:N_movie
    for j=1:2 % loop thoru channels
        % channel j
        xcor_img = zeros(img_size, img_size, n_rot);

        for r=1:n_rot % loop through rotations
            if j==1
                tmp = normxcorr2(lib_1(:,:,r), avg_img{i,j}); % x-correlate
            else
                tmp = normxcorr2(lib_2(:,:,r), avg_img{i,j}); % x-correlate
            end
            
            xcor_img(:,:,r) = tmp(box_size/2:end-box_size/2, box_size/2:end-box_size/2);
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
            p(:,1:2) = p(:,1:2)+box_size/2+1-1; % reshift and subrtract to obtian correct coordinates
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
%{
close all
particles = cell(N_movie, 2);
particles_rot = cell(N_movie, 2);
box_size_large = ceil(2*box_size/2/cos(pi/4));
dx = int16((box_size_large-box_size-1)/2);

for i=1:N_movie
    for j=2
    N = img_size;
    w = box_size/2;
    
    p = zeros(2*w, 2*w, size(peaks{i,j},1));
    p_rot = zeros(2*w, 2*w, size(peaks{i,j},1));

        for n=1:size(peaks{i,j},1) % loop through peak of image i,j
            x = peaks{i,j}(n,1);
            y = peaks{i,j}(n,2);
            angle = peaks{i,j}(n,4);
            p(max(1,w-y+2):min(2*w, 2*w+1-y-w+N-1),max(1,w-x+2):min(2*w, 2*w+1-x-w+N-1),j) = avg_img{i,j}(max(1, y-w):min(N,y+w-1) , max(1, x-w):min(N,x+w-1) );
            
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
            %imagesc(p_rot(:,:,j)), colorbar, colormap gray, axis image
            
            imagesc(lib_2(:,:,peaks{i,j}(n,5))), colorbar, colormap gray, axis image
            title(['Particle ' num2str(j) ' angle = ' num2str(angle) ' template ' num2str(peaks{i,j}(n,5))])
            pause
            
            
            
        end
        particles{i,j} = p;
        particles_rot{i,j} = p_rot;

    end
end
clear('p', 'p_rot')
%}
%%
%{
d = 4.8;
figure(2)
imagesc(avg_img{1,1}), axis image, colormap gray, hold on
plot(peaks{1,1}(:,1), peaks{1,1}(:,2), 'rx') 
dx = d*sin(peaks{1,1}(:,4)*pi/180)/2;
dy = d*cos(peaks{1,1}(:,4)*pi/180)/2;
plot(peaks{1,1}(:,1)-dx, peaks{1,1}(:,2)-dy, 'gx') 
plot(peaks{1,1}(:,1)+dx, peaks{1,1}(:,2)+dy, 'g+') 
%}
%% fiducial init_coordinates
d = 4.8;
coord_init = cell(N_movie, 2); 
for i=1:N_movie
    for j=1:2
        dx = d*sin(peaks{i,j}(:,4)*pi/180)/2;
        dy = d*cos(peaks{i,j}(:,4)*pi/180)/2;
        coord_init{i,j} = [peaks{i,j}(:,1) peaks{i,j}(:,2) peaks{i,j}(:,1)-dx peaks{i,j}(:,2)-dy peaks{i,j}(:,1)+dx peaks{i,j}(:,2)+dy ];
    end
end



%% fit fiducial positions
peaks_fit = cell(size(peaks));

for i=1:N_movie
    for j=1:2

        close all
        c_cur = zeros(size(peaks{i,j},1), 9);

        for n=1:size(coord_init{i,j},1)

            x0 = coord_init{i,j}(n,3);
            y0 = coord_init{i,j}(n,4);
            x1 = coord_init{i,j}(n,5);
            y1 = coord_init{i,j}(n,6);

            sigma_init = 1.5;
            w_fit = 8;
            [c, c_err, ci, area, residuals] = fit_2gauss2d_bg(x0, y0, x1, y1, sigma_init, w_fit, avg_img{i,j} ); %x1, y1, s_x, w_fit, images{peaks_raw(i, 5),1});

            c_cur(n, : ) = c;

            %{
            imagesc(avg_img{1,1}), axis image, colormap gray, hold on
              %  plot(coord_init{i,j}(n,1), coord_init{i,j}(n,2) , 'rx') 
            plot(x0, y0 , 'rx') 
            plot(x1, y1, 'r+')

            plot(c(1), c(2) , 'gx') 
            plot(c(5), c(6), 'g+') 

            [c(1:8) sqrt( (c(1)-c(5)).^2 + (c(2)-c(6)).^2 ) ]
            c_err(1:8)
            mse
            c(4)./c(8)
            pause

            hold off
            %}
        end
        peaks_fit{i,j} = c_cur;
    end
end


%% sort out single spots and bads spots
peaks_fit_filtered = cell(size(peaks));
peaks_filtered = cell(size(peaks));

plot_discarded = strcmp(questdlg('Plot discarded spots?','Plot discarded','Yes','No','No'), 'Yes');
if plot_discarded
    path_out_discarded = [path_out filesep 'discarded'];
    mkdir(path_out_discarded)
end

for i=1:N_movie
    for j=1:2
            N_peaks_raw = size(peaks{i,j},1);

        d_cutoff = 7;
        height_ratio = [0.2 2];

        criteria = ones(N_peaks_raw, 2);
        criteria(:,1) = sqrt((peaks_fit{i,j}(:,1)-peaks_fit{i,j}(:,5)).^2 + (peaks_fit{i,j}(:,2)-peaks_fit{i,j}(:,6)).^2) < d_cutoff;
        criteria(:,2) = height_ratio(1) < peaks_fit{i,j}(:,4)./peaks_fit{i,j}(:,8) & peaks_fit{i,j}(:,4)./peaks_fit{i,j}(:,8) < height_ratio(2);

        accepted = [criteria(:,1) & criteria(:,2)];

        %remove not-accepted spots
        peaks_filtered{i,j} = peaks{i,j}(accepted==1, :);
        peaks_fit_filtered{i,j} = peaks_fit{i,j}(accepted==1, :);




        close all
        fig_dim =1*[20 10];
        cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
        colormap gray
        w_plot = 10;


        if plot_discarded
            display('Discarding spots...')
            for k=1:N_peaks_raw
                if ~accepted(k) % discarded spot
                    %{
                    message = {'Sigma ratio green: OK', 'Sigma ratio red: OK', 'Spotsize green: OK', 'Spotsize red: OK'};
                    if criteria(i,1)==0
                        message{1} = 'Sigma ratio green: BAD';
                    end
                    if criteria(i,2)==0
                        message{2} = 'Sigma ratio red: BAD';
                    end
                    if criteria(i,3)==0
                        message{3} = 'Spotsize green: BAD';
                    end
                    if criteria(i,4)==0
                        message{4} = 'Spotsize red: BAD';
                    end
                    %}
                    x_1 = peaks_fit{i,j}(k,1);
                    y_1 = peaks_fit{i,j}(k,2);

                    x_2 = peaks_fit{i,j}(k,5);
                    y_2 = peaks_fit{i,j}(k,6);


                    plot_subframe(avg_img{i,j}, peaks{i,j}(k,1), peaks{i,j}(k,2), w_plot), hold on
                    plot(x_1, y_1, 'gx')
                    plot(x_2, y_2, 'g+')
                    ellipse(peaks_fit{i,j}(k,3), peaks_fit{i,j}(k,3), 0, x_1, y_1, 'g'  );
                    ellipse(peaks_fit{i,j}(k,7), peaks_fit{i,j}(k,7), 0, x_2, y_2, 'g'  );
                    %title({['Pair ' num2str(k) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in ' channel{1} ' channel'], message{1}, message{3}})
                    axis square
                    hold off

                    

                    print(cur_fig, '-dtiff', '-r150',  [path_out_discarded filesep 'Discarded_mov' num2str(i) '_ch' num2str(j) '_spot' num2str(k) '.tif'])
                end 
            end
        end

        display(['Disgarded ' num2str(sum(~accepted)) ' spot.'])

        close all

        N_peaks = size(peaks,1);
    end
end
    



%% find pairs
pairs = struct([]); 
for i=1:N_movie

    % map peaks
    trace_map = map_traces(peaks_filtered{i,1}(:, 1:2), peaks_filtered{i,2}(:, 1:2), peaks_filtered{i,1}(:, 1:2), 5)+1; %map the tarces from averga positions

    N_pairs = size(trace_map,1);

    pairs_tmp(1,N_pairs) = struct('red_coords', [], 'green_coords', [],  'red_cor', [],  'green_cor', [], 'movie', []);
    for k=1:size(trace_map,1)
        pairs_tmp(k).red_coords = peaks_fit_filtered{i,1}(trace_map(k,1),:); 
        pairs_tmp(k).green_coords = peaks_fit_filtered{i,2}(trace_map(k,2),:);
        pairs_tmp(k).red_cor = peaks_filtered{i,1}(trace_map(k,1),:);
        pairs_tmp(k).green_cor = peaks_filtered{i,2}(trace_map(k,2),:);
        pairs_tmp(k).movie = i;
    end
    pairs = [pairs, pairs_tmp];
    clear pairs_tmp
end

%% Mapping using fitgeotrans
tmp1 = [vertcat(pairs.red_coords) vertcat(pairs.movie)];
tmp2 = [vertcat(pairs.green_coords) vertcat(pairs.movie)];


xy_1 = [tmp1(:,1:2); tmp1(:,5:6) ];
xy_2 = [tmp2(:,1:2); tmp2(:,5:6)];
n_mov = [tmp1(:,10); tmp1(:,10)];

%tform_2TO1 = fitgeotrans(xy_2,xy_1,'polynomial',  2); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 1 = f(2)
%tform_1TO2 = fitgeotrans(xy_1,xy_2, 'polynomial',  2); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 2 = f(1)

tform_2TO1 = fitgeotrans(xy_2,xy_1,'lwm',  20); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 1 = f(2)
tform_1TO2 = fitgeotrans(xy_1,xy_2, 'lwm',  20); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 2 = f(1)


% Mapp coordinates
xy_tmp = transformPointsInverse(tform_2TO1, xy_1);  %%this takes coords in ch1 and transfroms them to coords in ch2
xy_ch2 = [xy_tmp xy_2 n_mov];

%sum(sum((xy_1-xy_2).^2))
%sum(sum(( transformPointsInverse(tform_2TO1, xy_1) -xy_2).^2))

xy_tmp = transformPointsInverse(tform_1TO2, xy_2);  %%this takes coords in ch2 and transfroms them to coords in ch1
xy_ch1 = [xy_1 xy_tmp n_mov];

xy_raw = [xy_1 xy_2 n_mov];
N_peaks = size(xy_ch1,1);


%%
N_pairs = size(pairs,2);
for i=1:N_pairs
    pairs(i).spot1_ch1 = [pairs(i).red_coords(1:2), transformPointsInverse(tform_1TO2, pairs(i).green_coords(1:2))];
    pairs(i).spot2_ch1 = [pairs(i).red_coords(5:6), transformPointsInverse(tform_1TO2, pairs(i).green_coords(5:6))];
    
    pairs(i).spot1_ch2 = [ transformPointsInverse(tform_2TO1, pairs(i).red_coords(1:2)), pairs(i).green_coords(1:2)];
    pairs(i).spot2_ch2 = [ transformPointsInverse(tform_2TO1, pairs(i).red_coords(5:6)) pairs(i).green_coords(5:6)];

end



%% SAVE everything
save([path_out filesep 'fiducial_mapping_data.mat'])
tform = tform_2TO1;
save([path_out filesep 'tform_ ' channel{2} '2' channel{1} '.mat'], 'tform')
tform = tform_1TO2;
save([path_out filesep 'tform_ ' channel{1} '2' channel{2} '.mat'], 'tform')
display(['Data saved to ' path_out])


%% Write images
path_out_img = [path_out filesep 'rgb_images'];
mkdir(path_out_img)

for i=1:N_movie
    rgb = zeros(size(avg_img{i,1},1), size(avg_img{i,1},1), 3, 'uint16');
        
    A = min_max_uint16(avg_img{i,1});
    B = min_max_uint16(avg_img{i,2});
    B_tformed= imwarp(B, tform_2TO1, 'OutputView', imref2d(size(A)));
    A_tformed= imwarp(A, tform_1TO2, 'OutputView', imref2d(size(B)));

   rgb(:,:,1) = A;
   rgb(:,:,2) = B_tformed;    
   imwrite(rgb, [path_out_img filesep 'image_' sprintf('%02i', i) '_on_' channel{1} '.tif'])
   
   
   rgb(:,:,1) = A_tformed;
   rgb(:,:,2) = B;    
   imwrite(rgb, [path_out_img filesep 'image_' sprintf('%02i', i) '_on_' channel{2} '.tif'] )
   
   rgb(:,:,1) = A;
   rgb(:,:,2) = B;    
   imwrite(rgb, [path_out_img filesep 'image_unmapped.tif'] )
end

% summed images
A = zeros(size(avg_img{1,1}));
B = zeros(size(avg_img{1,1}));
for i=1:N_movie     
    A = A + avg_img{i,1};
    B = B + avg_img{i,2};
end
A = min_max_uint16(A);
B = min_max_uint16(B);

B_tformed= imwarp(B, tform_2TO1, 'OutputView', imref2d(size(A)));
A_tformed= imwarp(A, tform_1TO2, 'OutputView', imref2d(size(B)));
rgb = zeros(size(avg_img{i,1},1), size(avg_img{i,1},1), 3, 'uint16');
rgb(:,:,1) = A;
rgb(:,:,2) = B_tformed;    
imwrite(rgb, [path_out_img filesep 'image_sum_' sprintf('%02i', i) '_on_' channel{1} '.tif'])


rgb(:,:,1) = A_tformed;
rgb(:,:,2) = B;    
imwrite(rgb, [path_out_img filesep 'image_sum_' sprintf('%02i', i) '_on_' channel{2} '.tif'] )

rgb(:,:,1) = A;
rgb(:,:,2) = B;    
imwrite(rgb, [path_out_img filesep 'image_sum_unmapped.tif'] )


%% plot each spot
close all
button = questdlg('Plot each pair?','Plot pair','Yes','No','No');
w_plot = 10;
if strcmp(button, 'Yes')
    path_out_pairs = [path_out filesep 'pairs'];
    mkdir(path_out_pairs)
    
    h = waitbar(0,'Printing spots... please wait');

    fig_dim =[20 10];
    cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    for i=1:N_pairs
 
                subplot(1, 2, 1) %image in channel 1
                x1 = pairs(i).spot1_ch1(1);
                y1 = pairs(i).spot1_ch1(2);
                
                x2 = pairs(i).spot2_ch1(1);
                y2 = pairs(i).spot2_ch1(2);
                
                plot_subframe(avg_img{pairs(i).movie, 1}, (x1+x2)/2, (y1+y2)/2, w_plot), hold on
                plot(x1, y1, 'rx', pairs(i).spot1_ch1(3), pairs(i).spot1_ch1(4), 'g+', 'MarkerSize', 15 )
                plot(x2, y2, 'rx', pairs(i).spot2_ch1(3), pairs(i).spot2_ch1(4), 'g+', 'MarkerSize', 15 )
                
                %plot(x, y, 'rx', xy_ch1(i,3), xy_ch1(i,4), 'g+', xy_raw(i,3), xy_raw(i,4), 'g.', 'MarkerSize', 15 )
                title(['Pair ' num2str(i) ' of '  num2str(N_pairs) ' at (' num2str(round(x1)) ',' num2str(round(y1)) ') in ' channel{1} ' channel'])
           %     legend({channel{1} , [channel{2} ' mapped'], [ channel{2} ' unmapped']})
                axis square
                hold off


                
                subplot(1, 2, 2) %image in channel 1
                x1 = pairs(i).spot1_ch2(1);
                y1 = pairs(i).spot1_ch2(2);
                
                x2 = pairs(i).spot2_ch2(1);
                y2 = pairs(i).spot2_ch2(2);
                
                plot_subframe(avg_img{pairs(i).movie, 2}, (x1+x2)/2, (y1+y2)/2, w_plot), hold on
                plot(x1, y1, 'g+', pairs(i).spot1_ch2(3), pairs(i).spot1_ch2(4), 'rx', 'MarkerSize', 15 )
                plot(x2, y2, 'g+', pairs(i).spot2_ch2(3), pairs(i).spot2_ch2(4), 'rx', 'MarkerSize', 15 )
                
                %plot(x, y, 'rx', xy_ch1(i,3), xy_ch1(i,4), 'g+', xy_raw(i,3), xy_raw(i,4), 'g.', 'MarkerSize', 15 )
                title(['Pair ' num2str(i) ' of '  num2str(N_pairs) ' at (' num2str(round(x1)) ',' num2str(round(y1)) ') in ' channel{2} ' channel'])
           %     legend({channel{1} , [channel{2} ' mapped'], [ channel{2} ' unmapped']})
                axis square
                hold off
                
               
                print(cur_fig, '-dtiff', '-r150',  [path_out_pairs filesep 'Pair_' num2str(i) '.tif'])
                
                waitbar( i/N_pairs , h, ['Printing spot... ' num2str(i) ' of ' num2str(N_pairs) ' done']) % update waitbar


    end             
    display('done printing pairs')
end
close(h)

%% finished
disp('DONE')




%%
i=1;
w_plot = 8;
close all
for k=1:length(pairs)
    subplot(2, 2, 1)
    plot_subframe(avg_img{i,1}, pairs(k).red_cor(1), pairs(k).red_cor(2), w_plot), axis square
    
    subplot(2, 2, 3)
    imagesc(lib_1(:,:,pairs(k).red_cor(5))), axis image

    subplot(2, 2, 2)
    plot_subframe(avg_img{i,2}, pairs(k).green_cor(1), pairs(k).green_cor(2), w_plot), axis square
    
    subplot(2, 2, 4)
    imagesc(lib_2(:,:,pairs(k).green_cor(5))), axis image

    pause
    
    
    
    
end



%%
for i=1:N_movie
    for j=1:2
    close all
    imagesc(avg_img{i,j}), axis image, colormap gray, hold on
    plot(peaks{i,j}(:,1), peaks{i,j}(:,2), 'rx') 
    plot(peaks_fit{i,j}(:,1), peaks_fit{i,j}(:,2), 'gx') 
    plot(peaks_fit{i,j}(:,5), peaks_fit{i,j}(:,6), 'g+') 
    title(['movie ' num2str(i) ' channel ' num2str(j)])

pause
    end
end
%%
%%
for i=1:N_movie
    for j=1:2
    close all
    imagesc(avg_img{i,j}), axis image, colormap gray, hold on
    plot(peaks_filtered{i,j}(:,1), peaks_filtered{i,j}(:,2), 'rx') 
    plot(peaks_fit_filtered{i,j}(:,1), peaks_fit_filtered{i,j}(:,2), 'gx') 
    plot(peaks_fit_filtered{i,j}(:,5), peaks_fit_filtered{i,j}(:,6), 'g+') 
    title(['movie ' num2str(i) ' channel ' num2str(j)])

pause
    end
end



