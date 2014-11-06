%% startup
clc, clear all, close all
path0 = cd;
run('my_prefs.m')


%% choose colors
rgb={'red','green','blue'};
[colors,ok]=listdlg('PromptString', 'Select ONE colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
while ne(length(colors),1)
    [colors,ok]=listdlg('PromptString', 'Select _TWO_ colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
end

channel = cell(1,1);
channel{1} = rgb{colors(1)};

%% LOAD STACK OF MOVIES
pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_ch1 = pickFirstFitsFiles(pname, channel{1}); 
N_movie = size(files_ch1,1);

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)

%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', ['Sequence ' channel{1} ':'],... % sample options
    'Radius of peak [pixel]:', 'Integration radius [pixel]:', 'Minimal length [frames]:'};
input_default = {'2', '-1', '1', '4', '3', '20'};
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

r_find = str2double(tmp(4)); % radius used to find spots
r_integrate = str2double(tmp(5)); % radius used for integration of intesituies
min_length = str2double(tmp(6)); % minimal number of found spots in a trace


%% generate movie classes
ch1 = cell(N_movie,1);

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
end

 

%% get threshold
% ch1 movie
for i=1:N_movie 
    ch1{i}.get_h_min(r_find);
end



%% trace movies
h = waitbar(0,'Traceing movies... please wait');
avg_img = cell(N_movie, 1);
positions = cell(N_movie, 1); % stores used spots
traces = cell(0,2);
itraces = cell(0,2);

for i=1:N_movie 
    [ch1_traces, ch1_itraces, avg_img{i,1}] = ch1{i}.trace_movie(ch1{i}.h_min, r_find, r_integrate, min_length);
    
    % create avergage position in ch1
    ch1_pos = zeros(size(ch1_traces,1),2);
    for j=1:size(ch1_traces,1)
        ch1_pos(j,1) = mean(ch1_traces{j,1}(:,2));
        ch1_pos(j,2) = mean(ch1_traces{j,1}(:,3));    
    end
    
    positions{i,1} = ch1_pos;    
    
    %get fluorescence traces from position
    ch1_itraces_full = ch1{i}.traces_movie_position(round(positions{i}(:,1:2)), r_integrate);
   
    movnumber = cell(size(ch1_itraces_full));
    movnumber(:) = {i};
    traces = [traces; ch1_traces movnumber]; % append traces to merged_traces
    itraces = [itraces; ch1_itraces_full movnumber];
    
    waitbar( i/N_movie , h, ['Traceing movies... ' num2str(i) ' of ' num2str(N_movie) ' done']) % update waitbar
end
close(h)


   
%% FIT GAUSSIAN TO DATA % MAP 1 On 2 and 2 ON 1 and store in merged_traces col 4:5
traces_fit = cell(size(traces));
h = waitbar(0,'Fitting PSF to spots... please wait');

for i=1:size(traces,1) % loop through traces   
    display(['Fitting spot ' num2str(i) ' of ' num2str(size(traces,1))])
    tmp = zeros(size(traces{i,1},1), size(traces{i,1},2)+9);
    tmp(:,1:3)= traces{i,1};
    tmp(:,6:12) = ch1{traces{i,2}}.fit_psf_to_movie(traces{i,1}(:,1), traces{i,1}(:,2:3), 2); %fit spot in each frame, sigma = 2
    traces_fit{i,1} = tmp;    
    
    waitbar( i/size(traces,1) , h, ['Fitting PSF to spots... ' num2str(i) ' of ' num2str(size(traces,1)) ' done']) % update waitbar

end
close(h)

%% save data
save([path_out filesep 'data.mat']);

%% Plot data
% traces for delta_r and sigma over frame number

path_out_traces = [path_out filesep 'position_traces'];
mkdir(path_out_traces)

cf = figure(1);
for i=1:size(traces,1)
    
    xy_c1 = traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    sigma_c1 = sqrt(traces_fit{i,1}(:,8).^2  + traces_fit{i,1}(:,9).^2);
    frame_c1 = traces_fit{i,1}(:,1);
    chi2_c1 = traces_fit{i,1}(:,12);
    meanchi2_c1=mean(chi2_c1);
    
    subplot(3, 1, 1)
    plot(frame_c1, d_c1, ['.-' rgb{colors(1)}(1)], 'Markersize', 15)
    legend(rgb{colors(1)})
    title(['Distance from average position. Spot: ' num2str(i)])
    xlabel('Frame')
    ylabel('$\sqrt{{\Delta x}^2 +{\Delta y}^2}$', 'Interpreter', 'latex')
    ylim([0 5])
    
    subplot(3, 1, 2)
    plot(frame_c1, sigma_c1, ['.-' rgb{colors(1)}(1)], 'Markersize', 15)
    legend(rgb{colors(1)})
    title(['Width of peak. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel('$\sqrt{{\sigma_x}^2+{\sigma_y}^2}$', 'Interpreter', 'latex')
    ylim([0 5])
    
    subplot(3, 1, 3)
    %bar(frame_green, chi2green)
    plot(frame_c1, chi2_c1, ['.-' rgb{colors(1)}(1)], 'Markersize', 15)
    title(['Chi^2 of ' rgb{colors(1)} ' fit. Mean value: ' sprintf('%.2s',meanchi2_c1) '. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel(['${\chi^2}_{' rgb{colors(1)} '}$'], 'Interpreter', 'latex')
    
    
    print(cf, '-dtiff','-r150', [path_out_traces filesep 'trace_' num2str(i) '.tif' ])
end

close all

%% plot intensity traces

path_itraces = [path_out filesep 'itraces'];
mkdir(path_itraces)

w_plot = 6;
close all
fig_dim =[20 10];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

h = waitbar(0,'Printing traces... please wait');

for i=1:size(itraces,1)
    I = itraces{i,1}(:,4);
    x = itraces{i,1}(1,2)+1;
    y = itraces{i,1}(1,3)+1;

    % image
    subplot(1,3,1)
    plot_subframe(avg_img{itraces{i,2},1}, x, y, w_plot), hold on
    plot(x, y, [channel{1}(1) '.'], 'MarkerSize', 15 )
    ellipse(r_integrate, r_integrate, 0, x, y, channel{1}(1));
    legend([num2str(x) ', ' num2str(y)])
    axis square
    hold off
    title(['Movie: ' num2str(itraces{i,2}) ', Trace: ' num2str(i) ' of ' num2str(size(itraces,1))])
 
    subplot(1,3, 2:3)
    plot(itraces{i,1}(:,1), itraces{i,1}(:,4), channel{1}(1),  'Linewidth', 1)
     ylabel('Integr. Intensity')
    
        
    print(cur_fig, '-dtiff', '-r150', [path_itraces filesep 'trace_' sprintf('%.03i',i) '.tif'])

    waitbar( i/size(itraces,1) , h, ['Printing traces... ' num2str(i) ' of ' num2str(size(itraces,1)) ' done']) % update waitbar

    
end
close(h)
close all

%% traces for position (dots unconnected)
close all
w_plot = 5;

cf = figure(1);
for i=1:size(traces,1)
    
    xy_c1 = traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    


    % image
    subplot(1,2,1)
    plot_subframe(avg_img{traces{i,2},1}, xy_mean_c1(1), xy_mean_c1(2), w_plot), hold on
    plot(xy_mean_c1(1), xy_mean_c1(2), [channel{1}(1) '.'], 'MarkerSize', 15 )
    legend([num2str(x) ', ' num2str(y)])
    axis square
    hold off
    title(['Movie: ' num2str(traces{i,2}) ', Trace: ' num2str(i) ' of ' num2str(size(traces,1))])
 
    
    subplot(1,2,2)
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) '.']), hold on
    xlim = int16([xy_mean_c1(1)-w_plot xy_mean_c1(1)+w_plot]);
    ylim = int16([xy_mean_c1(2)-w_plot xy_mean_c1(2)+w_plot]);
    set(gca, 'XLim', xlim, 'YLim', ylim, 'Ydir', 'reverse')
    axis square

    print(cf, '-dtiff', '-r150', [path_out_traces filesep 'position_' num2str(i) '.tif' ])
end

close all
%%


% traces for position (dots connected)
delta = 2;
cf = figure(1);
for i=1:size(traces,1)
    
    xy_c1 = traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) '.-']), hold on
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);

    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' fitted positions. Spot: ' num2str(i)]);
    hold off;
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
        axis square

    print(cf, '-dtiff', '-r150', [path_out_traces filesep 'position2_' num2str(i) '.tif' ])
end
hold off

close all
%%
% Plot positions of red and green spots on green average image
cf=figure(1);
for i=1:N_movie
    imagesc(avg_img{i,1}),colormap gray, colorbar;
    hold on
    plot(positions{i,1}(:,1),positions{i,1}(:,2),[channel{1}(1) 'o']);
    
    print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{1} '.eps'])
    hold off
end


%%
disp('Done')
