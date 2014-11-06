%% startup
clc, clear all, close all
path0 = cd;
run('my_prefs.m')


%% LOAD STACK OF MOVIES
pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_green = pickFirstFitsFiles(pname, 'green');  
files_red = pickFirstFitsFiles(pname, 'red');

N_movie = size(files_green,1);
if size(files_green,1) ~= size(files_red,1)
    disp('WARNING: not same number of movie files!')
end

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)
%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', 'Sequence green (D->D, D->A):', 'Sequence red (A->A):',... % sample options
    'Radius of peak [pixel]:', 'Integration radius [pixel]:', 'Minimal length [frames]:'};
input_default = {'2', '-1', '10', '01', '4', '3', '20'};
tmp = inputdlg(input, 'Parameters', 1, input_default);


first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file
%determine sequences 
sequence_dd = zeros(1, size(tmp{3},2));
sequence_da = zeros(1, size(tmp{3},2));
for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence_dd(1,i) =1;
        sequence_da(1,i) =1;
    end
end
sequence_aa = zeros(1, size(tmp{4},2));
for i=1:size(tmp{4},2)
    if(tmp{4}(i) == '1')
        sequence_aa(1,i) =1;
    end
end
r_find = str2double(tmp(5)); % radius used to find spots
r_integrate = str2double(tmp(6)); % radius used for integration of intesituies
min_length = str2double(tmp(7)); % minimal number of found spots in a trace



%% generate movie classes
dd = cell(N_movie,1);
da = cell(N_movie,1);
aa = cell(N_movie,1);

for i=1:N_movie
    dd{i} = movie(pname, files_green(i).name, first, last, sequence_dd); % pname, fname, first, last, sequence
    da{i} = movie(pname, files_red(i).name, first, last, sequence_da); % pname, fname, first, last, sequence
    aa{i} = movie(pname, files_red(i).name, first, last, sequence_aa); % pname, fname, first, last, sequence
end

%% get threshold
% dd movie
for i=1:N_movie 
    dd{i}.get_h_min(r_find);
end

%aa movie
for i=1:N_movie 
    aa{i}.get_h_min(r_find);
end

%% trace movies
h = waitbar(0,'Traceing movies... please wait');
avg_img = cell(N_movie, 3);
traces = cell(0, 4);
positions = cell(N_movie, 1);
mean_pos = zeros(0,6);

for i=1:N_movie 
    [dd_traces, dd_it, avg_img{i,1}] = dd{i}.trace_movie(dd{i}.h_min, r_find, r_integrate, min_length);
    
    %create avergage position of green
    dd_pos = zeros(size(dd_traces,1),2);
    for j=1:size(dd_traces,1)
        dd_pos(j,1) = mean(dd_traces{j,1}(:,2));
        dd_pos(j,2) = mean(dd_traces{j,1}(:,3));    
    end
    
    [aa_traces, aa_it, avg_img{i,3}] = aa{i}.trace_movie(aa{i}.h_min, r_find, r_integrate, min_length);
    %create avergage position of aa
    aa_pos = zeros(size(aa_traces,1),2);
    for j=1:size(aa_traces,1)
        aa_pos(j,1) = mean(aa_traces{j,1}(:,2));
        aa_pos(j,2) = mean(aa_traces{j,1}(:,3));    
    end
    
    avg_img{i,1} = dd{i}.average_image(dd{i}.mov_length);
    avg_img{i,2} = da{i}.average_image(da{i}.mov_length);
    avg_img{i,3} = aa{i}.average_image(aa{i}.mov_length);
    
    % map traces
    trace_map = map_traces(dd_pos, aa_pos, aa_pos, r_find)+1; %map the traces from averga positions
    positions{i} = zeros(size(trace_map,1), 6);
    for j=1:size(trace_map,1)
        positions{i}(j,:) = [dd_pos(trace_map(j,1),:) aa_pos(trace_map(j,2),:) aa_pos(trace_map(j,2),:)]; %%dd da aa
    end
    mean_pos = [mean_pos; positions{i}];
    
    %get fluorescence traces from position
    dd_it = dd{i}.traces_movie_position(positions{i}(:,1:2), r_integrate);
    da_it = da{i}.traces_movie_position(positions{i}(:,3:4), r_integrate);
    aa_it = aa{i}.traces_movie_position(positions{i}(:,5:6), r_integrate);
   
    movnumber = cell(size(dd_it));
    movnumber(:) = {i};
    traces = [traces; dd_it da_it aa_it movnumber];
    
    waitbar( i/N_movie , h, ['Traceing movies... ' num2str(i) ' of ' num2str(N_movie) ' done']) % update waitbar
end
close(h)


%% save data
save([path_out filesep 'data.mat'])
disp(['Saved data to ' path_out])

%%
d = sqrt(  (mean_pos(:,1)-mean_pos(:,3)).^2 + (mean_pos(:,2)-mean_pos(:,4)).^2  );



%{
w_plot = 6;
close all
for i=1:5%size(traces,1)
    Idd = traces{i,1}(:,4)-traces{i,1}(end,4);
    Ida = traces{i,2}(:,4)-traces{i,2}(end,4);
    
    subplot(2,3,1)
    x_dd = traces{i,1}(1,2)+1;
    y_dd = traces{i,1}(1,3)+1;
    
    plot_subframe(avg_img{traces{i,4},1}, x_dd, y_dd, w_plot), hold on
    plot(x_dd, y_dd, 'g.', 'MarkerSize', 15 )
    ellipse(r_integrate,r_integrate, 0, x_dd, y_dd, 'g')
    axis square
    hold off
    title(['Movie: ' num2str(traces{i,4}) ', Trace: ' num2str(i) ' of ' num2str(size(traces,1))])
     
    subplot(2,3, 2:3)
    plot(traces{i,1}(:,1), traces{i,1}(:,4)-traces{i,1}(end,4), 'g', traces{i,2}(:,1), traces{i,2}(:,4)-traces{i,2}(end,4), 'b', traces{i,3}(:,1), traces{i,3}(:,4)-traces{i,3}(end,4), 'r')
    
    subplot(2,3,4)
    x_aa = traces{i,3}(1,2)+1;
    y_aa = traces{i,3}(1,3)+1;
    
    plot_subframe(avg_img{traces{i,4},3}, x_aa, y_aa, w_plot), hold on
    plot(x_aa, y_aa, 'rx', 'MarkerSize', 15 )
    axis square
    hold off
     

    subplot(2,3,5:6)
    plot(traces{i,1}(:,1), Ida ./ (Idd+Ida), 'k')
    set(gca, 'YLim', [0 1])
    
    pause
end
%}

%% write txt output file

% number of frames in output matrix
N_min = size(traces{1,1},1); 
for i=1:size(traces,1)
    N_min = min([size(traces{i,1},1) size(traces{i,2},1) size(traces{i,3},1) N_min]);
end

for m=1:N_movie 
    floc= [path_out filesep 'data_movie_' sprintf('%.02i',m) '.txt']; % 

    index = find([traces{:,4}]==m); % index for specific movie
    N_tmp = length(index); % number of traces in this movie
    
    % generate output matrix
    A_out = zeros(N_min, 2+N_tmp*3);

    A_out(:,1) = traces{1,1}(1:N_min,1); % framenumber dd/da
    A_out(:,2) = traces{1,3}(1:N_min,1); % framenumber aa

    for i=1:N_tmp
        A_out(:,2+3*(i-1)+1) = traces{index(i),1}(1:N_min,4); %dd
        A_out(:,2+3*(i-1)+2) = traces{index(i),2}(1:N_min,4); %dd
        A_out(:,2+3*(i-1)+3) = traces{index(i),3}(1:N_min,4); %dd
    end

    % write output file
    % header file
    hdr = cell(1, 2+N_tmp*3);
    hdr{1} = 'frame_dd';
    hdr{2} = 'frame_aa';

     for i=1:N_tmp
         hdr{2+3*(i-1)+1} = ['dd_' num2str(index(i))];
         hdr{2+3*(i-1)+2} = ['da_' num2str(index(i))];
         hdr{2+3*(i-1)+3} = ['aa_' num2str(index(i))];
     end

    % write header file
    txt=sprintf('%s\t',hdr{:});
    txt(end)='';
    dlmwrite(floc,txt,'');
    dlmwrite(floc, A_out, '-append', 'delimiter', '\t'); % append data
end
disp('Data written to txt-file.')


%%
path_fig = [path_out filesep 'traces_fig'];
mkdir(path_fig)

w_plot = 6;
close all
fig_dim =[20 15];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

h = waitbar(0,'Printing traces... please wait');

for i=1:size(traces,1)
    Idd = traces{i,1}(:,4)-traces{i,1}(end,4);
    Ida = traces{i,2}(:,4)-traces{i,2}(end,4);
    
    x_dd = traces{i,1}(1,2)+1;
    y_dd = traces{i,1}(1,3)+1;
    x_da = traces{i,2}(1,2)+1;
    y_da = traces{i,2}(1,3)+1;
    x_aa = traces{i,3}(1,2)+1;
    y_aa = traces{i,3}(1,3)+1;
    
    
    % dd image
    subplot(3,3,1)
    plot_subframe(avg_img{traces{i,4},1}, x_dd, y_dd, w_plot), hold on
    plot(x_dd, y_dd, 'g.', 'MarkerSize', 15 )
    ellipse(r_integrate, r_integrate, 0, x_dd, y_dd, 'g');
    legend([num2str(x_dd) ', ' num2str(y_dd)])
    axis square
    hold off
    title(['Movie: ' num2str(traces{i,4}) ', Trace: ' num2str(i) ' of ' num2str(size(traces,1)) ', d = ' num2str(round(sqrt((x_dd-x_aa).^2 + (y_dd-y_aa).^2)*100)/100) 'px'])

    % da image
    subplot(3,3,2)
    plot_subframe(avg_img{traces{i,4},2}, x_da, y_da, w_plot), hold on
    plot(x_da, y_da, 'b.', 'MarkerSize', 15 )
    ellipse(r_integrate, r_integrate, 0 , x_da, y_da, 'b');
    legend([num2str(x_da) ', ' num2str(y_da)])
    axis square
    hold off
    
    % aa image
    subplot(3,3,3)
    plot_subframe(avg_img{traces{i,4},3}, x_aa, y_aa, w_plot), hold on
    plot(x_aa, y_aa, 'r.', 'MarkerSize', 15 )
    ellipse(r_integrate, r_integrate, 0, x_aa, y_aa, 'r');
    legend([num2str(x_aa) ', ' num2str(y_aa)])
    axis square
    hold off
    
    
    subplot(3,3, 4:6)
    plot(traces{i,1}(:,1), traces{i,1}(:,4), 'g', traces{i,2}(:,1), traces{i,2}(:,4), 'b', traces{i,3}(:,1), traces{i,3}(:,4), 'r', 'Linewidth', 1)
     ylabel('Integr. Intensity')
    
    
    subplot(3,3, 7:9)
    plot(traces{i,1}(:,1), Ida ./ (Idd+Ida), 'k', 'Linewidth', 1)
    xlabel('Framenumber'), ylabel('E_app')

    set(gca, 'YLim', [-0.1 1])
    
    print(cur_fig, '-dtiff', '-r300', [path_fig filesep 'trace_' sprintf('%.03i',i) '.tif'])

    waitbar( i/size(traces,1) , h, ['Printing traces... ' num2str(i) ' of ' num2str(size(traces,1)) ' done']) % update waitbar

    
end
close(h)
%% plot average images with traces
close all
fig_dim =2*[10 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

for i=1:N_movie
    
    
     %D_ex D_em
    imagesc(avg_img{i,1}),colormap gray, colorbar, axis image, hold on
    plot(positions{i}(:,1)+1 , positions{i}(:,2)+1  , 'go'), hold on
    title('D -> D'),
  
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image_' sprintf('%.02i',i) '_dd.tif'])
    
        
    
     %D_ex A_em
    imagesc(avg_img{i,2}),colormap gray, colorbar, axis image, hold on
    plot(positions{i}(:,1)+1 , positions{i}(:,2)+1  , 'bo'), hold on
    title('D -> A'),
  
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image_' sprintf('%.02i',i) '_da.tif'])
    
        
    
     %A_ex A_em
    imagesc(avg_img{i,3}),colormap gray, colorbar, axis image, hold on
    plot(positions{i}(:,1)+1 , positions{i}(:,2)+1  , 'ro'), hold on
    title('A -> A'),
  
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image_' sprintf('%.02i',i) '_aa.tif'])
end
%%
close all
display('Finished.')













