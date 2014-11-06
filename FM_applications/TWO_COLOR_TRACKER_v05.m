%% startup
clc, clear all, close all
path0 = cd;
run('my_prefs.m')


%% choose colors
rgb={'red','green','blue'};
[colors,ok]=listdlg('PromptString', 'Select two colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
while ne(length(colors),2)
    [colors,ok]=listdlg('PromptString', 'Select _TWO_ colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
end

channel = cell(2,1);
channel{1} = rgb{colors(1)};
channel{2} = rgb{colors(2)};

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
    'Radius of peak [pixel]:', 'Integration radius [pixel]:', 'Minimal length [frames]:'};
input_default = {'2', '-1', '1', '1', '4', '3', '20'};
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
r_integrate = str2double(tmp(6)); % radius used for integration of intesituies
min_length = str2double(tmp(7)); % minimal number of found spots in a trace


%% generate movie classes
ch1 = cell(N_movie,1);
ch2 = cell(N_movie,1);

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1(i).name, first, last, sequence_ch1); % pname, fname, first, last, sequence
    ch2{i} = movie(pname, files_ch2(i).name, first, last, sequence_ch2); % pname, fname, first, last, sequence
end

%%
button = questdlg(['Map positions ' channel{1} ' ON ' channel{2} ' and vice versa?'],'Mapping','Yes','No','No');
mapping = strcmp(button, 'Yes');
%%
if mapping == 1
    [mapping_file_1TO2, mapping_dir]=uigetfile(data_dir,['Choose the ' channel{1} '2' channel{2} ' mapping file:']);
    map1TO2 =load([mapping_dir mapping_file_1TO2], 'tform');
    tform_1TO2 = map1TO2.tform;
    display(['loaded ' channel{1} ' TO ' channel{2} ' mapping file: ' mapping_dir mapping_file_1TO2]);
    
    [mapping_file_2TO1]=uigetfile(mapping_dir,['Choose the ' channel{2} '2' channel{1} ' mapping file:']);
    map2TO1=load([mapping_dir mapping_file_2TO1]);
    tform_2TO1 = map2TO1.tform; %['tform_' channel{2} 'ON' channel{1}];
    display(['loaded ' channel{2} ' TO ' channel{1} ' mapping file: ' mapping_dir mapping_file_2TO1]);
end

 

%% get threshold
% ch1 movie
for i=1:N_movie 
    ch1{i}.get_h_min(r_find);
end

% ch2 movie
for i=1:N_movie 
    ch2{i}.get_h_min(r_find);
end



%% trace movies
h = waitbar(0,'Traceing movies... please wait');
avg_img = cell(N_movie, 2);
merged_traces = cell(0, 4);
positions = cell(N_movie, 1); % stores used spots
all_positions = cell(N_movie, 2); % stores all found spots
pos1on2 = cell(N_movie, 1);
pos2on1 = cell(N_movie, 1);

merged_traces = cell(0,3);
merged_itraces = cell(0,3);

for i=1:N_movie 
    [ch1_traces, ch1_itraces, avg_img{i,1}] = ch1{i}.trace_movie(ch1{i}.h_min, r_find, r_integrate, min_length);
    
    % create avergage position in ch1
    ch1_pos = zeros(size(ch1_traces,1),2);
    for j=1:size(ch1_traces,1)
        ch1_pos(j,1) = mean(ch1_traces{j,1}(:,2));
        ch1_pos(j,2) = mean(ch1_traces{j,1}(:,3));    
    end
    
    [ch2_traces, ch2_itraced, avg_img{i,2}] = ch2{i}.trace_movie(ch2{i}.h_min, r_find, r_integrate, min_length);
    %create avergage position in ch2
    ch2_pos = zeros(size(ch2_traces,1),2);
    for j=1:size(ch2_traces,1)
        ch2_pos(j,1) = mean(ch2_traces{j,1}(:,2));
        ch2_pos(j,2) = mean(ch2_traces{j,1}(:,3));    
    end
    
    if mapping
        pos1on2{i} = transformPointsInverse(tform_2TO1, ch1_pos);  %%this takes coords in ch1 and transfroms them to coords in ch2
        pos2on1{i} = transformPointsInverse(tform_1TO2, ch2_pos);  %%this takes coords in ch2 and transfroms them to coords in ch1
    end
    
    all_positions{i,1} = ch1_pos;
    all_positions{i,2} = ch2_pos;
    
    
    % map traces
    trace_map = map_traces(ch1_pos, ch2_pos, ch2_pos, r_find)+1; %map the traces from averga positions
    positions{i} = zeros(size(trace_map,1), 4);
    tmp = cell(size(trace_map,1),3);
    for j=1:size(trace_map,1)
        positions{i}(j,:) = [ch1_pos(trace_map(j,1),:) ch2_pos(trace_map(j,2),:) ]; %%dd da aa
        tmp{j,1} = ch1_traces{trace_map(j,1)};
        tmp{j,2} = ch2_traces{trace_map(j,2)};
        tmp{j,3} = i; % movie number
    end
    
    %get fluorescence traces from position
    ch1_itraces_full = ch1{i}.traces_movie_position(round(positions{i}(:,1:2)), r_integrate);
    ch2_itraces_full = ch2{i}.traces_movie_position(round(positions{i}(:,3:4)), r_integrate);
   
    movnumber = cell(size(ch1_itraces_full));
    movnumber(:) = {i};
    merged_traces = [merged_traces; tmp]; % append traces to merged_traces
    merged_itraces = [merged_itraces; ch1_itraces_full ch2_itraces_full movnumber];
    
    waitbar( i/N_movie , h, ['Traceing movies... ' num2str(i) ' of ' num2str(N_movie) ' done']) % update waitbar
end
close(h)


   
%% FIT GAUSSIAN TO DATA % MAP 1 On 2 and 2 ON 1 and store in merged_traces col 4:5
merged_traces_fit = cell(size(merged_traces));
h = waitbar(0,'Fitting PSF to spots... please wait');

for i=1:size(merged_traces,1) % loop through traces   
    display(['Fitting spot ' num2str(i) ' of ' num2str(size(merged_traces,1))])
    tmp = zeros(size(merged_traces{i,1},1), size(merged_traces{i,1},2)+9);
    tmp(:,1:3)= merged_traces{i,1};
    tmp(:,6:12) = ch1{merged_traces{i,3}}.fit_psf_to_movie(merged_traces{i,1}(:,1), merged_traces{i,1}(:,2:3), 2); %fit spot in each frame, sigma = 2
    if mapping == 1
        tmp(:,4:5) = transformPointsInverse(tform_2TO1, tmp(:,6:7));  %%this takes coords in ch1 and transfroms them to coords in ch2
    end
    merged_traces_fit{i,1} = tmp;    
    
    tmp = zeros(size(merged_traces{i,2},1), size(merged_traces{i,2},2)+9);
    tmp(:,1:3)= merged_traces{i,2};
    tmp(:,6:12) = ch2{merged_traces{i,3}}.fit_psf_to_movie(merged_traces{i,2}(:,1), merged_traces{i,2}(:,2:3), 2); %fit spot in each frame, sigma = 2
    if mapping == 1
        tmp(:,4:5) = transformPointsInverse(tform_1TO2, tmp(:,6:7));  %%this takes coords in ch2 and transfroms them to coords in ch1
    end
    merged_traces_fit{i,2} = tmp;   
    
    waitbar( i/size(merged_traces,1) , h, ['Fitting PSF to spots... ' num2str(i) ' of ' num2str(size(merged_traces,1)) ' done']) % update waitbar

end
close(h)

%% save data
save([path_out filesep 'data.mat']);

%% Plot data
% traces for delta_r and sigma over frame number
cf = figure(1);
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    sigma_c1 = sqrt(merged_traces_fit{i,1}(:,8).^2  + merged_traces_fit{i,1}(:,9).^2);
    frame_c1 = merged_traces_fit{i,1}(:,1);
    chi2_c1 = merged_traces_fit{i,1}(:,12);
    meanchi2_c1=mean(chi2_c1);
    
    xy_c2 = merged_traces_fit{i,2}(:,6:7);
    xy_mean_c2 = [mean(xy_c2(:,1)) mean(xy_c2(:,2))  ];
    d_c2 = sqrt((xy_c2(:,1)-xy_mean_c2(1)).^2 + (xy_c2(:,2)-xy_mean_c2(2)).^2);
    
    sigma_c2 = sqrt(merged_traces_fit{i,2}(:,8).^2  + merged_traces_fit{i,2}(:,9).^2);
    frame_c2 = merged_traces_fit{i,2}(:,1);
    chi2_c2 = merged_traces_fit{i,2}(:,12);
    meanchi2_c2 = mean(chi2_c2);

    subplot(4, 1, 1)
    plot(frame_c1, d_c1, ['.-' rgb{colors(1)}(1)],frame_c2,d_c2, ['.-' rgb{colors(2)}(1)], 'Markersize', 15)
    legend(rgb{colors(1)}, rgb{colors(2)})
    title(['Distance from average position. Spot: ' num2str(i)])
    xlabel('Frame')
    ylabel('$\sqrt{{\Delta x}^2 +{\Delta y}^2}$', 'Interpreter', 'latex')
    ylim([0 5])
    
    subplot(4, 1, 2)
    plot(frame_c1, sigma_c1, ['.-' rgb{colors(1)}(1)],frame_c2,sigma_c2, ['.-' rgb{colors(2)}(1)], 'Markersize', 15)
    legend(rgb{colors(1)}, rgb{colors(2)})
    title(['Width of peak. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel('$\sqrt{{\sigma_x}^2+{\sigma_y}^2}$', 'Interpreter', 'latex')
    ylim([0 5])
    
    subplot(4, 1, 3)
    %bar(frame_green, chi2green)
    plot(frame_c1, chi2_c1, ['.-' rgb{colors(1)}(1)], 'Markersize', 15)
    title(['Chi^2 of ' rgb{colors(1)} ' fit. Mean value: ' sprintf('%.2s',meanchi2_c1) '. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel(['${\chi^2}_{' rgb{colors(1)} '}$'], 'Interpreter', 'latex')
    
    subplot(4, 1, 4)
    plot(frame_c2, chi2_c2, ['.-' rgb{colors(2)}(1)], 'Markersize', 15)
    title(['Chi^2 of ' rgb{colors(2)} ' fit. Mean value: ' sprintf('%.2s',meanchi2_c2) '. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel(['${\chi^2}_{' rgb{colors(2)} '}$'], 'Interpreter', 'latex')
    
    print(cf, '-depsc2', [path_out filesep 'trace_' num2str(i) '.eps' ])
end

close all

%% traces for position (dots unconnected)
delta = 5;
cf = figure(1);
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    xy_c2 = merged_traces_fit{i,2}(:,6:7);
    xy_mean_c2 = [mean(xy_c2(:,1)) mean(xy_c2(:,2))  ];
    d_c2 = sqrt((xy_c2(:,1)-xy_mean_c2(1)).^2 + (xy_c2(:,2)-xy_mean_c2(2)).^2);
    
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) 'o']), hold on
    plot(xy_c2(:,1), xy_c2(:,2), [rgb{colors(2)}(1) 'x']), hold off
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' and ' rgb{colors(2)} ' fitted positions (unmapped). Spot: ' num2str(i)]);
    hold off;
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    print(cf, '-depsc2', [path_out filesep 'position_' num2str(i) '.eps' ])
end

close all


%% traces for position (dots unconnected, red dots mapped rONg)
if mapping == 1
delta = 5;
cf = figure(1);
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    xy_2on1 = merged_traces_fit{i,2}(:,4:5);
    xy_mean_2on1 = [mean(xy_2on1(:,1)) mean(xy_2on1(:,2))  ];
    d_rONg = sqrt((xy_2on1(:,1)-xy_mean_2on1(1)).^2 + (xy_2on1(:,2)-xy_mean_2on1(2)).^2);
    
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) 'o']), hold on
    plot(xy_2on1(:,1), xy_2on1(:,2), [rgb{colors(2)}(1) 'x']), hold off
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' and ' rgb{colors(2)} ' fitted positions (rONg mapped). Spot: ' num2str(i)]);
    hold off;
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    print(cf, '-depsc2', [path_out filesep 'position_mapped_rONg_' num2str(i) '.eps' ])
end
end

close all

%% traces for position (dots connected)
delta = 5;
cf = figure(1);
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    xy_c2 = merged_traces_fit{i,2}(:,6:7);
    xy_mean_c2 = [mean(xy_c2(:,1)) mean(xy_c2(:,2))  ];
    d_c2 = sqrt((xy_c2(:,1)-xy_mean_c2(1)).^2 + (xy_c2(:,2)-xy_mean_c2(2)).^2);
    
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) '-o']), hold on
    plot(xy_c2(:,1), xy_c2(:,2), [rgb{colors(2)}(1) '-x']), hold off
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' and ' rgb{colors(2)} ' fitted positions (unmapped). Spot: ' num2str(i)]);
    hold off;
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    print(cf, '-depsc2', [path_out filesep 'position2_' num2str(i) '.eps' ])
end
close all
hold off

%%
% Plot positions of red and green spots on green average image
cf=figure(1);
for i=1:N_movie
    imagesc(avg_img{i,1}),colormap gray, colorbar;
    hold on
    plot(all_positions{i,1}(:,1),all_positions{i,1}(:,2),[channel{1}(1) 'o']);
    plot(all_positions{i,2}(:,1),all_positions{i,2}(:,2),[channel{2}(1) 'x']);
    
    print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{1} '_' channel{2} '.eps'])
    hold off
end

% Plot positions of red and green spots on green average image
cf=figure(1);
for i=1:N_movie
    imagesc(avg_img{i,2}),colormap gray, colorbar;
    hold on
    plot(all_positions{i,1}(:,1),all_positions{i,1}(:,2),[channel{1}(1) 'o']);
    plot(all_positions{i,2}(:,1),all_positions{i,2}(:,2),[channel{2}(1) 'x']);
    
    print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{2} '_' channel{1} '.eps'])
    hold off
end


%%
if mapping == 1
    
    cf=figure(1);
    for i=1:N_movie
        imagesc(avg_img{i,1}),colormap gray, colorbar;
        hold on
        plot(all_positions{i,1}(:,1),all_positions{i,1}(:,2),[channel{1} 'o']);
        plot(pos2on1{i}(:,1),pos2on1{i}(:,2),[channel{2} 'x']);

        print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{2} 'TO' channel{1} '.eps'])
        hold off
    end
    
    for i=1:N_movie
        imagesc(avg_img{i,2}),colormap gray, colorbar;
        hold on
        plot(pos1on2{i}(:,1),pos1on2{i}(:,2),[channel{1} 'o']);
        plot(all_positions{i,2}(:,1),all_positions{i,2}(:,2),[channel{2} 'x']);

        print(cf, '-depsc2', [path_out filesep 'mov_' sprintf('%02i', i) '_spots_' channel{1} 'TO' channel{2} '.eps'])
        hold off
    end
    
end


%%
disp('Done')
