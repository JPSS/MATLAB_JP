%% STARTUP
clc, clear all, close all
path0 = cd; 
run('my_prefs')

%% Select channels
h.f = figure('units','pixels','position',[scrsz(3)/2,scrsz(4)/2,230,90],...
             'toolbar','none','menu','none');
% Create title
h.t = uicontrol('style','text','units','pixels','position',[10,70,210,15],'string','Select 2 channels and hit ENTER');
% Create yes/no checkboxes
h.c(1) = uicontrol('style','checkbox','units','pixels',...
                'position',[50,45,50,15],'string','red');
h.c(2) = uicontrol('style','checkbox','units','pixels',...
                'position',[50,30,50,15],'string','green');
h.c(3) = uicontrol('style','checkbox','units','pixels',...
                'position',[50,15,50,15],'string','blue');
pause;
vars = get(h.c,'Value');
close all

channel = cell(2,1);
j=1;
if vars{1}
    channel{j} = 'red'; j=j+1; 
end
if vars{2}
    channel{j} = 'green'; j=j+1;
end
if vars{3}
    channel{j} = 'blue'; j=j+1;
end
%% LOAD STACK OF MOVIES/FILES
pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_ch1 = pickFirstFitsFiles(pname, channel{1});
files_ch2 = pickFirstFitsFiles(pname, channel{2});

N_movie = size(files_ch1,1); % number of movies
if size(files_ch1,1) ~= size(files_ch2, 1)
    disp('WARNING: not same number of movie files!')
end

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_mapping'];
mkdir(path_out)
%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', ['Sequence ' channel{1} ':'], ['Sequence ' channel{2} ':'],... % sample options
    'Radius of peak [pixel]:', 'Integration radius [pixel]:', 'Frames for average (-1=all):'};
input_default = {'2', '-1', '1', '1', '4', '3', '-1'};
tmp = inputdlg(input, 'Parameters', 1, input_default);

first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file
%determine sequences 
sequence_1 = zeros(1, size(tmp{3},2));
for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence_1(1,i) =1;
    end
end
sequence_2 = zeros(1, size(tmp{4},2));
for i=1:size(tmp{4},2)
    if(tmp{4}(i) == '1')
        sequence_2(1,i) =1;
    end
end
r_find = str2double(tmp(5)); % radius used to find spots
r_integrate = str2double(tmp(6)); % radius used for integration of intesituies
N_frames = str2double(tmp(7)); % radius used for integration of intesituies

%% generate movie classes
ch1 = cell(N_movie,1);
ch2 = cell(N_movie,1);  

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1(i).name, first, last, sequence_1); % pname, fname, first, last, sequence
    ch2{i} = movie(pname, files_ch2(i).name, first, last, sequence_2); % pname, fname, first, last, sequence
end

%% determine thresholds and find peaks
peaks_raw = zeros(0,5);

for i=1:N_movie
    [h_min, p1] = ch1{i}.get_h_min(r_find, N_frames);
    [h_min, p2] = ch2{i}.get_h_min(r_find, N_frames);
    
    % map peaks
    trace_map = map_traces(p1(:,1:2), p2(:,1:2), p2(:,1:2), r_find*2)+1; %map the tarces from averga positions

    tmp = zeros(size(trace_map,1),5);
    % combine pairs
    for j=1:size(trace_map,1)
        tmp(j,:) = [p1(trace_map(j,1), 1:2)+1 p2(trace_map(j,2), 1:2)+1 i]; %x_1 y_1 x_2 y_2 frame
    end
    
    peaks_raw = [peaks_raw; tmp];
end

N_peaks_raw = size(peaks_raw,1);
display(['You have ' num2str(N_peaks_raw) ' pairs'])

%% compute averages images
images = cell(N_movie, 2);

for i=1:N_movie
    images{i, 1} = ch1{i}.average_image(N_frames);
    images{i, 2} = ch2{i}.average_image(N_frames);
end

%% Fit psf to spots
s_x = 2.5;
s_y = 2.5;
w_fit = 10;

ch1_fit_raw = zeros(N_peaks_raw, 7); 
ch1_fit_err_raw = zeros(N_peaks_raw, 7); 
ch2_fit_raw = zeros(N_peaks_raw, 7); 
ch2_fit_err_raw = zeros(N_peaks_raw, 7); 

h = waitbar(0,'Fitting spots... please wait');

for i=1:N_peaks_raw
    %display(['Fitting spot ' num2str(i) ' of ' num2str(size(peaks,1))])
    
    % channel 1
    x1 = round(peaks_raw(i,1));
    y1 = round(peaks_raw(i,2));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x1, y1, s_x, w_fit, images{peaks_raw(i, 5),1});
    ch1_fit_raw(i,:) = c;
    ch1_fit_err_raw(i,:) = c_err;
    
    % channel 2
    x2 = round(peaks_raw(i,3));
    y2 = round(peaks_raw(i,4));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x2, y2, s_x, w_fit, images{peaks_raw(i, 5),2});
    ch2_fit_raw(i,:) = c;
    ch2_fit_err_raw(i,:) = c_err;
    
    waitbar( i/N_peaks_raw , h, ['Fitting spot... ' num2str(i) ' of ' num2str(N_peaks_raw) ' done']) % update waitbar

end

close(h)

%% SORT OUT: remove spots where ratio of width is not close to 1 and which are to large

criteria = ones(N_peaks_raw,4 );
criteria(:,1:2) = filter_spots(ch1_fit_raw(:,3:4), [0.8 1.2], 2);
criteria(:,3:4) = filter_spots(ch2_fit_raw(:,3:4), [0.8 1.2], 2);
accepted = [criteria(:,1) & criteria(:,2) & criteria(:,3) & criteria(:,4)];

%remove not-accepted spots
ch1_fit = ch1_fit_raw(accepted==1, :);
ch1_fit_err = ch1_fit_err_raw(accepted==1, :);
ch2_fit = ch2_fit_raw(accepted==1, :);
ch2_fit_err = ch2_fit_err_raw(accepted==1, :);
peaks = peaks_raw(accepted==1, :);

plot_discarded = strcmp(questdlg('Plot discarded spots?','Plot discarded','Yes','No','No'), 'Yes');
if plot_discarded
    path_out_discarded = [path_out filesep 'discarded'];
    mkdir(path_out_discarded)
end

close all
fig_dim =1*[20 10];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
colormap gray
w_plot = 10;


if plot_discarded
    display('Discarding spots...')
    for i=1:N_peaks_raw
        if ~accepted(i) % discarded spot
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
            x_1 = ch1_fit_raw(i,1);
            y_1 = ch1_fit_raw(i,2);

            x_2 = ch2_fit_raw(i,1);
            y_2 = ch2_fit_raw(i,2);

            subplot(1, 2, 1)
            plot_subframe(images{peaks_raw(i, 5), 1}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'rx')
            plot(x_2, y_2, 'g+')
            ellipse(ch1_fit_raw(i,3), ch1_fit_raw(i,4), -ch1_fit_raw(i,5), x_1, y_1, 'r'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in ' channel{1} ' channel'], message{1}, message{3}})
            axis square
            hold off

            subplot(1, 2, 2)
            plot_subframe(images{peaks_raw(i, 5), 2}, x_2, y_2, w_plot), hold on
            plot(x_1, y_1, 'rx')
            plot(x_2, y_2, 'g+')
            ellipse(ch2_fit_raw(i,3), ch2_fit_raw(i,4), -ch2_fit_raw(i,5), x_2, y_2, 'r'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_2)) ',' num2str(round(y_2)) ') in ' channel{2} ' channel'], message{2}, message{4}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_discarded filesep 'Discarded_' num2str(i) '.tif'])
        end 
    end
end

display(['Disgarded ' num2str(sum(~accepted)) ' spot.'])

close all

N_peaks = size(peaks,1);
%% Mapping using fitgeotrans
xy_1 = ch1_fit(:,1:2);
xy_2 = ch2_fit(:,1:2);

tform_2TO1 = fitgeotrans(xy_2,xy_1,'lwm', 6); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 1 = f(2)
tform_1TO2 = fitgeotrans(xy_1,xy_2,'lwm', 6); %moving_points, fixed_points, local weihted mean, nearest naeighbours, 2 = f(1)

%% Mapp coordinates
xy_tmp = transformPointsInverse(tform_2TO1, xy_1);  %%this takes coords in ch1 and transfroms them to coords in ch2
xy_ch2 = [xy_tmp xy_2 peaks(:,5)];

%sum(sum((xy_1-xy_2).^2))
%sum(sum(( transformPointsInverse(tform_2TO1, xy_1) -xy_2).^2))

%%
xy_tmp = transformPointsInverse(tform_1TO2, xy_2);  %%this takes coords in ch2 and transfroms them to coords in ch1
xy_ch1 = [xy_1 xy_tmp peaks(:,5)];

xy_raw = [xy_1 xy_2 peaks(:,5)];

%% SAVE everything
save([path_out filesep 'data.mat'])
tform = tform_2TO1;
save([path_out filesep 'tform_ ' channel{2} '2' channel{1} '.mat'], 'tform')
tform = tform_1TO2;
save([path_out filesep 'tform_ ' channel{1} '2' channel{2} '.mat'], 'tform')
display(['Data saved to ' path_out])

%% Write images
path_out_img = [path_out filesep 'rgb_images'];
mkdir(path_out_img)

for i=1:N_movie
    rgb = zeros(size(images{i,1},1), size(images{i,1},1), 3, 'uint16');
        
    A = min_max_uint16(images{i,1});
    B = min_max_uint16(images{i,2});
    B_tformed= imwarp(B, tform_2TO1, 'OutputView', imref2d(size(A)));
    A_tformed= imwarp(A, tform_1TO2, 'OutputView', imref2d(size(B)));

   rgb(:,:,1) = A;
   rgb(:,:,2) = B_tformed;    
   imwrite(rgb, [path_out_img filesep 'image_' sprintf('%02i', i) '_on_' channel{1} '.tif'])
   
   
   rgb(:,:,1) = A_tformed;
   rgb(:,:,2) = B;    
   imwrite(rgb, [path_out_img filesep 'image_' sprintf('%02i', i) '_on_' channel{2} '.tif'] )
end


%% PLOT distances
display('------------------------- PLOTTING DATA --------------------------')

% plot distance before and after mapping
index = 1:N_peaks;
d = sqrt((xy_raw(:,1)-xy_raw(:,3)).^2 + (xy_raw(:,2)-xy_raw(:,4)).^2);
d_tform = sqrt(   (xy_ch1(:,1)-xy_ch1(:,3)).^2 + (xy_ch1(:,2)-xy_ch1(:,4)).^2   );
chi2_tform = sum(d_tform.^2);

close all
fig_dim =2*[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(1,4,[1:3])
plot(index, d, '.-k', index, d_tform, '.-r','MarkerSize', 15)
legend({'unmapped', 'mapped'})
xlabel('Spot')
ylabel('Distance')

subplot(1,4,4)
bar([1 2], [sum(d.^2) sum(d_tform.^2) ])
set(gca, 'XTick', [1:2] , 'XTickLabel', {'Unmapped', 'mapped'})
set(gca, 'XLim', [0 3])
ylabel('chi2')

print(cur_fig, '-dtiff', '-r300', [path_out filesep 'Mapping_comparison.tif'])


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
    for i=1:N_peaks
 
                subplot(1, 2, 1) %image in channel 1
                x = xy_ch1(i,1);
                y = xy_ch1(i,2);
                plot_subframe(images{xy_ch1(i,5), 1}, x, y, w_plot), hold on
                plot(x, y, 'rx', xy_ch1(i,3), xy_ch1(i,4), 'g+', xy_raw(i,3), xy_raw(i,4), 'g.', 'MarkerSize', 15 )
                title(['Pair ' num2str(i) ' of '  num2str(N_peaks) ' at (' num2str(round(x)) ',' num2str(round(y)) ') in ' channel{1} ' channel'])
                legend({channel{1} , [channel{2} ' mapped'], [ channel{2} ' unmapped']})
                axis square
                hold off

                subplot(1, 2, 2) %image in channel 2
                x = xy_ch2(i,3);
                y = xy_ch2(i,4);
                plot_subframe(images{xy_ch2(i,5), 2}, x, y, w_plot), hold on
                plot(x, y, 'gx', xy_ch2(i,1), xy_ch2(i,2), 'r+', xy_raw(i,1), xy_raw(i,2), 'r.', 'MarkerSize', 15 )
                title(['Pair ' num2str(i) ' of '  num2str(N_peaks) ' at (' num2str(round(x)) ',' num2str(round(y)) ') in ' channel{2} ' channel'])
                legend({channel{2} , [channel{1} ' mapped'], [ channel{1} ' unmapped']})
                axis square
                hold off

               
                print(cur_fig, '-dtiff', '-r150',  [path_out_pairs filesep 'Pair_' num2str(i) '.tif'])
                
                waitbar( i/N_peaks , h, ['Printing spot... ' num2str(i) ' of ' num2str(N_peaks) ' done']) % update waitbar


    end             
    display('done printing pairs')
end
close(h)



%% print distance scatterplot and get width
close all
fig_dim =[20 20];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

dx = xy_ch1(:,3)-xy_ch1(:, 1);
dy = xy_ch1(:,4)-xy_ch1(:, 2);
lim = 1.1*[min([min(dx) min(dy)]) max([max(dx) max(dy)]) ];

xhist = lim(1):(lim(2)-lim(1))/10:lim(2);

n_dx = hist(dx, xhist);
n_dy = hist(dy, xhist);


subplot(2, 2, 1)
plot(dx, dy, 'r.', 'MarkerSize', 5)
set(gca, 'XLim', lim)
set(gca, 'YLim', lim)
axis square
ylabel('Distance [pixel]')

subplot(2, 2, 3)
bar(xhist, n_dx, 'r'), hold on
set(gca, 'XLim', lim)
set(gca, 'Ydir', 'reverse')
axis square
legend(['sigma_x= ' num2str(std(dx)) ])
xlabel('Distance [pixel]')

subplot(2, 2, 2)
barh(xhist, n_dy, 'r'), hold on
legend(['sigma_y= ' num2str(std(dy)) ])
set(gca, 'YLim', lim)
axis square

print(cur_fig, '-dtiff', '-r300', [path_out filesep 'Distance_scatter_' channel{1} '.tif'])






%% distance distribution
close all

d = sqrt((xy_raw(:,1)-xy_raw(:,3)).^2 + (xy_raw(:,2)-xy_raw(:,4)).^2);
d_tform = sqrt(   (xy_ch1(:,1)-xy_ch1(:,3)).^2 + (xy_ch1(:,2)-xy_ch1(:,4)).^2   );


xhist = 0:max(d)/20:max(d);
n = hist(d, xhist);

xhist2 = 0:max(d)/200:max(d);
n_tform = hist(d_tform, xhist2);

fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(2, 1, 1)
bar(xhist, n)
set(gca, 'XLim', [xhist(1) xhist(end)])
ylabel('Frequency')
legend({'no correction'})
title('Channel 1')


subplot(2, 1, 2)
bar(xhist2, n_tform)
set(gca, 'XLim', [xhist(1) xhist(end)])
xlabel('Distance [pixel]')
ylabel('Frequency')
legend({'mapped'})
print(cur_fig, '-dtiff', '-r300', [path_out filesep 'Distance_distribution_' channel{1} '.tif'])


%% plot histogram of width sigma 

fig_dim =[15 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

xhist = 0:0.05:4;
[n_ch1] = hist(sqrt(ch1_fit(:,3).^2 + ch1_fit(:,4).^2),xhist);
[n_ch2] = hist(sqrt(ch2_fit(:,3).^2 + ch2_fit(:,4).^2),xhist);
subplot(3,1, 1)
plot(xhist, n_ch1, '-r', xhist, n_ch2, '-g', 'MarkerSize', 15)
set(gca, 'XLim', [xhist(1) xhist(end)])
title('Width of peak')
ylabel('Frequency')
legend({[channel{1} ' channel'], [ channel{2} ' channel']})

subplot(3,1, 2)
bar(xhist, n_ch1, 'r')
set(gca, 'XLim', [xhist(1) xhist(end)])
ylabel('Frequency')
legend({[channel{1} ' channel']})

subplot(3, 1, 3)
bar(xhist, n_ch2, 'g')
set(gca, 'XLim', [xhist(1) xhist(end)])
ylabel('Frequency')
legend({[channel{2} ' channel']})

ylabel('Frequency')
xlabel('$\sqrt{\sigma_x^2 + \sigma_y^2}$ in pixel', 'Interpreter', 'LaTex')

print(cur_fig, '-dtiff', '-r300',  [path_out filesep 'width.tif'])



%% scatterplot of width
close all
fig_dim =[10 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
plot(ch1_fit(:,3), ch1_fit(:,4), 'r.', ch2_fit(:,3), ch2_fit(:,4), 'g.', 'MarkerSize', 15)
xlabel('$\sigma_x$ in pixel', 'Interpreter', 'LaTex')
ylabel('$\sigma_y$ in pixel', 'Interpreter', 'LaTex')
title('Width - scatterplot')
legend({[channel{1} ' channel'], [channel{2} ' channel']})

lim = [0.9*min(min([ch1_fit(:,3) ch1_fit(:,4) ch2_fit(:,3) ch2_fit(:,4)])) 1.1*max(max([ch1_fit(:,3) ch1_fit(:,4) ch2_fit(:,3) ch2_fit(:,4)]))];


set(gca, 'Xlim', lim)
set(gca, 'Ylim', lim)
axis square
print(cur_fig, '-dtiff', '-r300',   [path_out filesep 'width_scatterplot.tif'])




%% print a vector-plot , unscaled
fig_dim =[30 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(1, 2, 1)
dx = xy_ch2(:,1)-xy_ch1(:, 1);
dy = xy_ch2(:,2)-xy_ch1(:, 2);
quiver(xy_ch1(:,1), xy_ch1(:,2), dx, dy, channel{1}(1) ,'AutoScale','off')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{1} ' positions'])

subplot(1, 2, 2)
dx = xy_ch1(:,3)-xy_ch2(:, 3);
dy = xy_ch1(:,4)-xy_ch2(:, 4);
quiver(xy_ch2(:,3), xy_ch2(:,4), dx, dy, channel{2}(1) ,'AutoScale','off')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{2} ' positions'])

print(cur_fig, '-depsc2', [path_out filesep 'Transformation_Image_absolute.eps'])

%% print a vector-plot , scaled
fig_dim =[30 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(1, 2, 1)
dx = xy_ch2(:,1)-xy_ch1(:, 1);
dy = xy_ch2(:,2)-xy_ch1(:, 2);
quiver(xy_ch1(:,1), xy_ch1(:,2), dx, dy, channel{1}(1) ,'AutoScale','on')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{1} ' positions'])

subplot(1, 2, 2)
dx = xy_ch1(:,3)-xy_ch2(:, 3);
dy = xy_ch1(:,4)-xy_ch2(:, 4);
quiver(xy_ch2(:,3), xy_ch2(:,4), dx, dy, channel{2}(1) ,'AutoScale','on')
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{2} ' positions'])

print(cur_fig, '-depsc2', [path_out filesep 'Transformation_Image_scaled.eps'])


%%
[X,Y] = meshgrid(50:10:450,50:10:450);
xy = [X(:) Y(:)];


close all
fig_dim =[30 30];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

xy_map = transformPointsInverse(tform_2TO1, xy);  %%this takes coords in ch1 and transfroms them to coords in ch2

dx = xy_map(:,1)-xy(:, 1);
dy = xy_map(:,2)-xy(:, 2);
quiver(xy(:,1), xy(:,2), dx, dy, channel{1}(1) ,'AutoScale','off', 'LineWidth', 1)
set(gca, 'XLim', [-20 size(images{1,1},1)+20])
set(gca, 'YLim', [-20 size(images{1,1},1)+20])
axis image
title(['Transformation of ' channel{1} ' to ' channel{2}])

print(cur_fig, '-depsc2', [path_out filesep 'Transformation_Image_grid.eps'])

%%
close all
display('Mappping finished.')
