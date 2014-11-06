%% startup
clc, clear all, close all
path0 = cd;
run('my_prefs.m')


%% LOAD STACK OF MOVIES
pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_green = pickFirstFitsFiles(pname, 'green');
files_red = pickFirstFitsFiles(pname, 'red');

N_movie = size(files_green,2);
if size(files_green,1) ~= size(files_red,1)
    disp('WARNING: not same number of movie files!')
end

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)
%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', 'Sequence green (D->D, D->A):', 'Sequence red (A->A):',... % sample options
    'Radius of peak [pixel]:', 'Integration radius [pixel]:', 'Average Frame length [frames]:'};
input_default = {'2', '-1', '10', '01', '4', '3', '-1'};
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
N_frames = str2double(tmp(7)); % minimal number of found spots in a trace



%% generate movie classes
dd = cell(N_movie,1);
da = cell(N_movie,1);
aa = cell(N_movie,1);

for i=1:N_movie
    dd{i} = movie(pname, files_green{i}, first, last, sequence_dd); % pname, fname, first, last, sequence
    da{i} = movie(pname, files_red{i}, first, last, sequence_da); % pname, fname, first, last, sequence
    aa{i} = movie(pname, files_red{i}, first, last, sequence_aa); % pname, fname, first, last, sequence
end


%% determine thresholds and find peaks
peaks_raw = zeros(0,5);

for i=1:N_movie
    [h_min, pdd] = dd{i}.get_h_min(r_find, N_frames);
    [h_min, paa] = aa{i}.get_h_min(r_find, N_frames);
    
    % map peaks
    trace_map = map_traces(pdd(:,1:2), paa(:,1:2), paa(:,1:2), r_find*2)+1; %map the tarces from averga positions

    tmp = zeros(size(trace_map,1),5);
    % combine pairs
    for j=1:size(trace_map,1)
        tmp(j,:) = [pdd(trace_map(j,1), 1:2)+1 paa(trace_map(j,2), 1:2)+1 i]; %x_1 y_1 x_2 y_2 frame
    end
    
    peaks_raw = [peaks_raw; tmp];
end

N_peaks_raw = size(peaks_raw,1);
display(['You have ' num2str(N_peaks_raw) ' pairs'])

%% compute averages images
avg_img = cell(N_movie, 3);

for i=1:N_movie
    avg_img{i, 1} = dd{i}.average_image(N_frames);
    avg_img{i, 2} = da{i}.average_image(N_frames);
    avg_img{i, 3} = aa{i}.average_image(N_frames);
end


%% Fit psf to spots
s_x = 2.5;
s_y = 2.5;
w_fit = 10;

dd_fit_raw = zeros(N_peaks_raw, 7); 
dd_fit_err_raw = zeros(N_peaks_raw, 7); 
aa_fit_raw = zeros(N_peaks_raw, 7); 
aa_fit_err_raw = zeros(N_peaks_raw, 7); 

h = waitbar(0,'Fitting spots... please wait');

for i=1:N_peaks_raw
    %display(['Fitting spot ' num2str(i) ' of ' num2str(size(peaks,1))])
    
    % channel 1
    x1 = round(peaks_raw(i,1));
    y1 = round(peaks_raw(i,2));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x1, y1, s_x, w_fit, avg_img{peaks_raw(i, 5),1});
    dd_fit_raw(i,:) = c;
    dd_fit_err_raw(i,:) = c_err;
    
    % channel 2
    x2 = round(peaks_raw(i,3));
    y2 = round(peaks_raw(i,4));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x2, y2, s_x, w_fit, avg_img{peaks_raw(i, 5),3});
    aa_fit_raw(i,:) = c;
    aa_fit_err_raw(i,:) = c_err;
    
    waitbar( i/N_peaks_raw , h, ['Fitting spot... ' num2str(i) ' of ' num2str(N_peaks_raw) ' done']) % update waitbar

end

close(h)

%% SORT OUT: remove spots where ratio of width is not close to 1 and which are to large

criteria = ones(N_peaks_raw,4 );
criteria(:,1:2) = filter_spots(dd_fit_raw(:,3:4), [0.8 1.2], 2);
criteria(:,3:4) = filter_spots(aa_fit_raw(:,3:4), [0.8 1.2], 2);
accepted = [criteria(:,1) & criteria(:,2) & criteria(:,3) & criteria(:,4)];

%remove not-accepted spots
dd_fit = dd_fit_raw(accepted==1, :);
dd_fit_err = dd_fit_err_raw(accepted==1, :);
aa_fit = aa_fit_raw(accepted==1, :);
aa_fit_err = aa_fit_err_raw(accepted==1, :);
peaks = peaks_raw(accepted==1, :);
peaks = [dd_fit(:,1:2) aa_fit(:,1:2) peaks(:,5)]; % use fitted poistions for further analysis


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
                message{1} = 'Sigma ratio D->D: BAD';
            end
            if criteria(i,2)==0
                message{2} = 'Sigma ratio A->A: BAD';
            end
            if criteria(i,3)==0
                message{3} = 'Spotsize D->D: BAD';
            end
            if criteria(i,4)==0
                message{4} = 'Spotsize A->A: BAD';
            end
            x_1 = dd_fit_raw(i,1);
            y_1 = dd_fit_raw(i,2);

            x_2 = aa_fit_raw(i,1);
            y_2 = aa_fit_raw(i,2);

            subplot(1, 2, 1)
            plot_subframe(avg_img{peaks_raw(i, 5), 1}, x_1, y_1, w_plot), hold on
            plot(x_1, y_1, 'g.')
            ellipse(dd_fit_raw(i,3), dd_fit_raw(i,4), -dd_fit_raw(i,5), x_1, y_1, 'g'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_1)) ',' num2str(round(y_1)) ') in D->D channel'], message{1}, message{3}})
            axis square
            hold off

            subplot(1, 2, 2)
            plot_subframe(avg_img{peaks_raw(i, 5), 3}, x_1, y_1, w_plot), hold on
            plot(x_2, y_2, 'r.')
            ellipse(aa_fit_raw(i,3), aa_fit_raw(i,4), -aa_fit_raw(i,5), x_2, y_2, 'r'  );
            title({['Pair ' num2str(i) ' of '  num2str(N_peaks_raw) ' at (' num2str(round(x_2)) ',' num2str(round(y_2)) ') in A->A channel'], message{2}, message{4}})
            axis square
            hold off

            print(cur_fig, '-dtiff', '-r150',  [path_out_discarded filesep 'Discarded_' num2str(i) '.tif'])
        end 
    end
end

display(['Disgarded ' num2str(sum(~accepted)) ' spot.'])
close all
N_peaks = size(peaks,1);


%% trace movies
h = waitbar(0,'Traceing movies... please wait');
traces = cell(0, 4);
positions = cell(N_movie, 1);

for i=1:N_movie 
    
    dd_pos = peaks(peaks(:,5)==i, 1:2); % peaks for this movie
    aa_pos = peaks(peaks(:,5)==i, 3:4); % peaks for this movie

    positions{i} = [dd_pos aa_pos i*ones(size(dd_pos,1),1)];
    
    %get fluorescence traces from position
    dd_it = dd{i}.traces_movie_position(round(dd_pos(:,1:2)-1), r_integrate);
    da_it = da{i}.traces_movie_position(round(aa_pos(:,1:2)-1), r_integrate);
    aa_it = aa{i}.traces_movie_position(round(aa_pos(:,1:2)-1), r_integrate);
   
    movnumber = cell(size(dd_it));
    movnumber(:) = {i};
    traces = [traces; dd_it da_it aa_it movnumber];
    
    waitbar( i/N_movie , h, ['Traceing movies... ' num2str(i) ' of ' num2str(N_movie) ' done']) % update waitbar
end
close(h)


%% save data
save([path_out filesep 'data.mat'])
disp(['Saved data to ' path_out])

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

d = sqrt(  (peaks(:,1)-peaks(:,3)).^2 + (peaks(:,2)-peaks(:,4)).^2  );

path_fig = [path_out filesep 'traces_fig'];
mkdir(path_fig)

w_plot = 6;
close all
fig_dim =[20 15];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

h = waitbar(0,'Printing traces... please wait');

for i=1:size(traces,1)
    %Idd = traces{i,1}(:,4)-traces{i,1}(end,4);
    %Ida = traces{i,2}(:,4)-traces{i,2}(end,4);
    
    Idd = traces{i,1}(:,4);
    Ida = traces{i,2}(:,4);
    
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
    plot(traces{i,3}(:,1), traces{i,3}(:,4)./max(traces{i,3}(:,4)), 'r', traces{i,1}(:,1), traces{i,1}(:,4)./max(traces{i,1}(:,4)), 'g', traces{i,2}(:,1), traces{i,2}(:,4)./max(traces{i,2}(:,4)), 'b',  'Linewidth', 1)
     ylabel('Integr. Intensity')
    
    
    subplot(3,3, 7:9)
    plot(traces{i,1}(:,1), Ida ./ (Idd+Ida), 'k', 'Linewidth', 1)

    xlabel('Framenumber'), ylabel('E_app')
    set(gca, 'YLim', [min(Ida ./ (Idd+Ida)) max(Ida ./ (Idd+Ida))])
    %set(gca, 'YLim', [-0.1 1])
    
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
    plot(positions{i}(:,1) , positions{i}(:,2)  , 'go'), hold on
    title('D -> D'),
  
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image_' sprintf('%.02i',i) '_dd.tif'])
    
        
    
     %D_ex A_em
    imagesc(avg_img{i,2}),colormap gray, colorbar, axis image, hold on
    plot(positions{i}(:,3) , positions{i}(:,4)  , 'bo'), hold on
    title('D -> A'),
  
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image_' sprintf('%.02i',i) '_da.tif'])
    
        
    
     %A_ex A_em
    imagesc(avg_img{i,3}),colormap gray, colorbar, axis image, hold on
    plot(positions{i}(:,3) , positions{i}(:,4)  , 'ro'), hold on
    title('A -> A'),
  
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image_' sprintf('%.02i',i) '_aa.tif'])
end
%%
close all
display('Finished.')













