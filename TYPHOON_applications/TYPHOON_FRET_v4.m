%% startup
clc, clear all, close all
path0 = cd; addpath(path0); display(['Added search-path: ' path0 ])
run('my_prefs')

%% Load images
cd(data_dir)
[filename_dd pathname_dd]=uigetfile('*.tif','Select Dex -> Dem [cy3]image: ');
cd(pathname_dd)
[filename_da pathname_da]=uigetfile('*.tif','Select Dex -> Aem [FRET] image: ');
[filename_aa pathname_aa]=uigetfile('*.tif','Select Aex -> Aem [cy5] image: ');
cd(path0)

dd_bin = double(imread([pathname_dd filesep filename_dd])); %bin4x4(double(imread([pathname_dd filesep filename_dd])));
da_bin = double(imread([pathname_da filesep filename_da])); %bin4x4(double(imread([pathname_da filesep filename_da])));
aa_bin = double(imread([pathname_aa filesep filename_aa])); %bin4x4(double(imread([pathname_aa filesep filename_aa])));

%% check PMT values of images
pmt = zeros(3,1);
dd_info = imfinfo([pathname_dd filesep filename_dd]); 
i_found = strfind(dd_info.ImageDescription, 'PMT=');
disp([ 'D->D channel [cy3] PMT: ' dd_info.ImageDescription(i_found+9:i_found+9+3)])
pmt(1) = str2double(dd_info.ImageDescription(i_found+9:i_found+9+2));

da_info = imfinfo([pathname_da filesep filename_da]); 
i_found = strfind(da_info.ImageDescription, 'PMT=');
disp([ 'D->A channel [FRET] PMT: ' da_info.ImageDescription(i_found+9:i_found+9+3)])
pmt(2) = str2double(da_info.ImageDescription(i_found+9:i_found+9+2));

aa_info = imfinfo([pathname_aa filesep filename_aa]); 
i_found = strfind(aa_info.ImageDescription, 'PMT=');
disp([ 'A->A channel [cy5] PMT: ' aa_info.ImageDescription(i_found+9:i_found+9+3)])
pmt(3) = str2double(aa_info.ImageDescription(i_found+9:i_found+9+2));

%% check for saturation
%generate a uint16 colormap
b = [0:1:(2^16-1)]/(2^16-1);
cm_uint16 = [b' b' b'];
cm_uint16(end,:) = [1 0 0 ];

close all
if max(dd_bin(:)) == 2^16-1
    disp('WARNING: D->D image saturated')
    imagesc(dd_bin), colormap(cm_uint16), colorbar, axis image, hold on
    [x, y] = find(dd_bin == 2^16-1);
    plot(y, x, 'r.')
    questdlg('WARNING: D->D image saturated','Saturation','Ignore','Ignore');
else
    disp('D->D image good')
end    
if max(da_bin(:)) == 2^16-1
    disp('WARNING: D->A image saturated')
    imagesc(da_bin), colormap(cm_uint16), colorbar, axis image, hold on
    [x, y] = find(da_bin == 2^16-1);
    plot(y, x, 'r.')
    questdlg('WARNING: D->A image saturated','Saturation','Ignore','Ignore');
else
    disp('D->A image good')
end
if max(aa_bin(:)) == 2^16-1
    disp('WARNING: A->A image saturated')
    imagesc(aa_bin), colormap(cm_uint16), colorbar, axis image, hold on
    [x, y] = find(aa_bin == 2^16-1);
    plot(y, x, 'r.')
    questdlg('WARNING: A->A image saturated','Saturation','Ignore','Ignore');
else
    disp('A->A image good')
end
close all

%% create output folder
pname = inputdlg({'Output folder and prefix:'}, 'Output folder and prefix' , 1, {[filename_dd(1:end-10) '_analysis']} );
prefix_out = pname{1};
path_out = [pathname_dd prefix_out ];
mkdir(path_out);
path_out_plots = [path_out filesep 'plots'];
mkdir(path_out_plots)

%%  correct for background
[dd_bg, bg_dd ] = bg_correct_ui(dd_bin, 'Donor excitation -> Donor emission');
[da_bg,bg_da ] = bg_correct_ui(da_bin, 'Donor excitation -> Acceptor emission');
[aa_bg, bg_aa ] = bg_correct_ui(aa_bin, 'Acceptor excitation -> Acceptor emission');

%% shift images
[dd_bg, da_x_min, da_y_min] = overlay_image(da_bg, dd_bg, 10);
[aa_bg, aa_x_min, aa_y_min]= overlay_image(da_bg, aa_bg, 10);

%% leakage and direct-excitation correction factors
close all
button = questdlg('Use leak/dir-corrections from this gel or load old data?','Leak/Dir','This Gel','Load Data','No correction', 'Load Data');
correction = 1;
if strcmp(button,'This Gel') %load old data
    leak_dir  = calculate_corrections(dd_bg, da_bg, aa_bg, [path_out filesep prefix_out '_correction.txt']);
    da_cor = da_bg - leak_dir(1,1).*dd_bg - leak_dir(2,1).*aa_bg;%-leak_dir(1,2)-leak_dir(2,2); 
    display('Used this gel to correct.')
else
    if strcmp(button,'Load Data')
        cd(data_dir)
        [fname_cor pname_cor] = uigetfile('.txt', 'Choose file with correction matrix');
        leak_dir = load([pname_cor fname_cor]);
        cd(path0)
        da_cor = da_bg - leak_dir(1,1).*dd_bg - leak_dir(2,1).*aa_bg; 
        display('Used old data to correct.')
    else
        display('No correction will be done.')
        correction = 0;
        da_cor = da_bg; 
        leak_dir = [1 1 ; 1 1 ];
    end
end


%% plot uncorrected and corrected images
scrsz = get(0,'ScreenSize');
cur_fig = figure('Visible','on','OuterPosition',[ 1 scrsz(4) scrsz(3) scrsz(4)/1.5], 'PaperPositionMode', 'auto'); % figure('Visible','off');%left bottom width height
subplot(1,2,1)
imagesc(da_bg), axis image, colormap gray
title('D->A image, Uncorrected')
subplot(1,2,2)
imagesc(da_cor), axis image, colormap gray
title('D->A image, Leak/Dir-corrected')
pause(2)
close all

%% gamma correction
button = questdlg('Do gamma-correction?','Gamma correction','This Gel','Load gamma','No','No');
if strcmp(button,'This Gel') 
    % selct reference for gamma-factor determination
    cd(pathname_dd)
    [filename_ref, pathname_ref]=uigetfile('*.tif','Select reference image (Cy2): ');
    cd(path0)
    ref_raw = -double(imread([pathname_ref filesep filename_ref])); 
    ref_bg = bg_correct_ui(ref_raw, 'Reference image');


    gamma = estimate_gamma_bands(dd_bg, da_cor, aa_bg, ref_bg) ; 
    save([path_out filesep prefix_out '_gamma.txt'], 'gamma' ,'-ascii')


else
    if strcmp(button,'Load gamma') 
        cd(data_dir)
        [fname_gamma pname_gamma] = uigetfile('.txt', 'Choose file with gamma-factor');
        gamma = load([pname_gamma fname_gamma]);
        cd(path0)     
    else
        gamma =1 ;
    end
end

%% find lanes
[auto_pos , area] = find_lanes(da_bg+aa_bg+dd_bg);
area = [min(area(:,1)) min(area(:,2)) max(area(:,1)+area(:,3))-min(area(:,1)) max(area(:,2)+area(:,4))-min(area(:,2))]; %area which surrounds all subareas
n_lanes = size(auto_pos,1);
%% generate profiles lanes

lanes = cell(0, 7);
for i=1:size(auto_pos,1)

    new_lane = cell(1, 7);

    pos = auto_pos(i, :);
    y = transpose(pos(2):pos(2)+pos(4));
    new_lane{1, 1} = +1 .* transpose( sum(transpose(  dd_bg(y  ,   pos(1):pos(1)+pos(3)   )  )) );
    new_lane{1, 2} =    +1 .* transpose( sum(transpose(  da_cor(y  ,   pos(1):pos(1)+pos(3)   )  )) );
    new_lane{1, 3} =    +1 .* transpose( sum(transpose(  aa_bg(y  ,   pos(1):pos(1)+pos(3)   )  )) );
    new_lane{1, 4} = y;
    new_lane{1,5} = pos;

    new_lane{1,6} = ['Lane ' num2str(size(lanes,1)+1)]; 
    new_lane{1,7} = i+9; 
    lanes = [lanes; new_lane];

end
%% Naming lanes
close all
button = questdlg('Name lanes?','Name lanes','Yes','No','No');
if strcmp(button,'Yes') %load old data
    for i=1:n_lanes
        
        subplot(2, 1, 1)
        imagesc(da_cor), axis image, colormap gray
        title('D->A')
        hold on
        rectangle('Position', auto_pos(i,:), 'EdgeColor', 'r')
        
        subplot(2, 1, 2)
        imagesc(da_cor), axis image, colormap gray
        title('D->A')
        hold on
        for j=1:i-1
            rectangle('Position', lanes{j,5}, 'EdgeColor', 'r')
            text(double(lanes{j,5}(1)), double(lanes{j,5}(2)), lanes{j,6}, 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'bottom')
        end
        hold off

        lane_name = inputdlg({'Name of lane:', 'Length of spacer:'}, 'Lane properties' , 1, {['Lane ' num2str(1)], num2str(i+9)} );
        lanes{1,7} = str2double(lane_name{2});
        lanes{1,6} = lane_name{1};
        
        display([ lane_name{1} ':   spacer=' lane_name{2}])
    end
end


%% fit gaussian close to maxima
options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',
bands = cell(0,5); % pos, label, spacer lenbth
bands_fit = cell(size(lanes,1),3);

I_sum = zeros(size(lanes,1), 3);
I_max = zeros(size(lanes,1), 3);
y_max = zeros(size(lanes,1), 3);
fit_dd = zeros(n_lanes, 4); fit_da = zeros(n_lanes, 4); fit_aa = zeros(n_lanes, 4);
sum_limits = zeros(n_lanes, 2);

close all
for i=1:size(lanes,1)

    % Fit a gaussian to the profile
    % display(['Fitting ' lanes{i,6}])
    p_dd = fit_lane(double(lanes{i,4}), lanes{i,1}, 5, 1, 1, 0); % pixel, intensity, initial width of peak, sigma_left, sigma_right, no bg
    p_da = fit_lane(double(lanes{i,4}), lanes{i,2}, 5, 1, 1, 0); % pixel, intensity, initial width of peak, sigma_left, sigma_right, no bg
    p_aa = fit_lane(double(lanes{i,4}), lanes{i,3}, 5, 1, 1, 0); % pixel, intensity, initial width of peak, sigma_left, sigma_right, no bg
    
    % save fitted parameters
    fit_dd(i,:) = p_dd;
    fit_da(i,:) = p_da;
    fit_aa(i,:) = p_aa;
   
    % filter da signal for maximum determination
    sigma_filter = 2;
    data = lanes{i,2};
    threeSigma = ceil(3*sigma_filter);
    g = exp(-(-threeSigma:threeSigma).^2/2/sigma_filter^2);
    da = conv(data, g, 'same') ./ conv(ones(size(data)), g, 'same');
    
    % FIND PEAKS IN DA-channel
    peaks_da = find_peaks1d(da, 20, 0.5*max(da), 1); % window size 20, h_min=0.5*max, 1 = absolute height
    y_mean = peaks_da(end);  % use the last found peak (leading band)
    dy = 5; % integration width
    
    % sum lanes
    I_sum(i,1) = sum(lanes{i,1}(y_mean-dy:y_mean+dy));
    I_sum(i,2) = sum(lanes{i,2}(y_mean-dy:y_mean+dy));
    I_sum(i,3) = sum(lanes{i,3}(y_mean-dy:y_mean+dy));
    
    sum_limits(i,:) = [y_mean dy];

    I_max(i,1) = p_dd(3); % height + bg
    I_max(i,2) = p_da(3);
    I_max(i,3) = p_aa(3);
    
    y_max(i,1) = p_dd(1);  
    y_max(i,2) = p_da(1);
    y_max(i,3) = p_aa(1);


    pos =[ lanes{i,5}(1) y_mean-dy+lanes{i,5}(2)-1 lanes{i,5}(3) 2*dy];
    
    bands{i,1} = -1*[ sum(sum(dd_bg( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)))) sum(sum(da_cor( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)))) sum(sum(aa_bg( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)))) ; sum(sum(dd_bg( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)))) sum(sum(da_bg( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)))) sum(sum(aa_bg( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3))));   sum(sum(dd_bin( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)))) sum(sum(da_bin( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3) ))) sum(sum(aa_bin( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)))) ];
    bands{i,2} = pos ; 
    bands{i,3} = lanes{i,6};
    bands{i,4} = lanes{i,7};
    

    
   %plot(lanes{i,4},lanes{i,1},'g',  lanes{i,4},lanes{i,2},'b',  lanes{i,4},lanes{i,3},'r' ), hold on
   %plot(y, da, 'b--')
   %plot(y, gauss1d(p_dd, y), 'g--', y, gauss1d(p_da, y), 'b--', y, gauss1d(p_aa, y), 'r--')
   %vline(y(y_mean-dy), 'k--')
  % vline(y(y_mean+dy), 'k--')
  % vline(y(y_mean), 'k')
 %  hold off
   %pause
   %pause(1)
end


close all

%% write corrected images
disp('Writing images')
close all
%imwrite(uint16(da_cor), [path_out filesep 'da_cor.tif'] , 'tif')
%imwrite(uint16(da_cor-min(min(da_cor))), [path_out filesep 'da_cor2.tif'] , 'tif')
%imwrite(uint16(-da_cor+abs(min(min(-da_cor)))), [path_out filesep 'da_cor3.tif'] , 'tif')

t = Tiff([path_out filesep 'da_cor.tif'],'w');
t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
t.setTag('BitsPerSample',16);
t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
t.setTag('ImageLength',size(da_cor,1));
t.setTag('ImageWidth',size(da_cor,2));
t.setTag('SamplesPerPixel',1);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write( uint16(da_cor-min(da_cor(:)))  );
t.close();

t = Tiff([path_out filesep 'da_cor_2.tif'],'w');
t.setTag('Photometric',Tiff.Photometric.MinIsBlack);
t.setTag('BitsPerSample',16);
t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
t.setTag('ImageLength',size(da_cor,1));
t.setTag('ImageWidth',size(da_cor,2));
t.setTag('SamplesPerPixel',1);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write( uint16(da_cor-min(da_cor(:)))  );
t.close();

t = Tiff([path_out filesep 'da_cor+bg.tif'],'w');
t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
t.setTag('BitsPerSample',16);
t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
t.setTag('ImageLength',size(da_cor,1));
t.setTag('ImageWidth',size(da_cor,2));
t.setTag('SamplesPerPixel',1);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write( uint16(da_cor+bg_da)  );
t.close();



%% save all the data
disp('Saving data...')
save([path_out filesep prefix_out '_data.mat'])

%% sava data for fret only
save([path_out filesep prefix_out '_data_FRET.mat'], 'I_sum', 'gamma')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generating plots...')
%% Plot profiles seperately 
path_out_profiles = [path_out filesep 'profiles'];
mkdir(path_out_profiles)

%determine y limits
ylim = [ 0 0 ];
for i=1:size(lanes,1)
    if max(max([lanes{i,1} lanes{i,2} lanes{i,3}])) > ylim(2)
        ylim(2) = max(max([lanes{i,1} lanes{i,2} lanes{i,3}]));
    end
    if min(min([lanes{i,1} lanes{i,2} lanes{i,3}])) > ylim(1)
        ylim(1) = min(min([lanes{i,1} lanes{i,2} lanes{i,3}]));
    end
    
end
ylim(2) = 1.1*ylim(2);

close all    
fig_dim =[20 8];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

h = zeros(3, 1);
for i=1:size(lanes,1)

    h(1:3) = plot(lanes{i,4},lanes{i,1},'g',  lanes{i,4},lanes{i,2},'b',  lanes{i,4},lanes{i,3},'r' ); hold on
   % h(4:6) = plot(lanes{i,4}, gauss1d( fit_dd(i,:) , lanes{i,4}), 'g--' , lanes{i,4}, gauss1d( fit_da(i,:) , lanes{i,4}), 'b--'  , lanes{i,4}, gauss1d( fit_aa(i,:) , lanes{i,4}), 'r--'   );
    xlabel('Pixel')
   ylabel('Intensity [a.u.]')
      
   title([lanes{i,6} ' (' num2str(i) ' of ' num2str(size(lanes, 1)) ')'])
   set(gca, 'XLim', [min(lanes{i,4}) max(lanes{i,4})])
   set(gca, 'YLim', ylim)
   vline(lanes{i,4}(sum_limits(i,1) +  sum_limits(i,2)), 'k--');
   vline(lanes{i,4}(sum_limits(i,1) -  sum_limits(i,2)), 'k--');
   %legend(h, {'D -> D', 'D -> A', 'A -> A', 'D -> D, fit', 'D -> A, fit', 'A -> A, fit'})
   legend(h, {'D -> D', 'D -> A', 'A -> A'})
   print(cur_fig, '-dtiff', '-r300', [path_out_profiles filesep 'profile_lane_' sprintf('%.02i',i) '.tif'])

   hold off
end
close all


%% plot FRET-profile of all lanes
close all
fig_dim =[20 8];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
for i=1:size(lanes,1)
    E = lanes{i,2} ./ (lanes{i,2} + gamma * lanes{i,1});
   for j=1:size(E,1)
       if lanes{i,2}(j) < 0
           E(j) = 0;
       end
       
       if lanes{i,1}(j) < 0
           E(j) = 0;
       end
   end
   h(i) = plot(lanes{i,4}, E, 'Linewidth', 1);
end
set(gca, 'YLim', [0 1])
legend(h, lanes{:,6} )
xlabel('Pixel')
ylabel('FRET efficiency')
title({'FRET efficiency', ['gamma = ' num2str(gamma)]})
print(cur_fig, '-dtiff' , '-r500', [path_out_plots filesep 'FRET_profiles.tif']); %save figure

  
%% plot migration distance and width of peak

close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

width = zeros(size(lanes,1), 3);
y_max_max = zeros(size(lanes,1), 3);

for i=1:size(lanes,1)
    width(i,:) = [fit_dd(i, 2) fit_da(i, 2) fit_aa(i, 2) ];
    [tmp j1] =max(lanes{i,1});
    [tmp j2] =max(lanes{i,2});
    [tmp j3] =max(lanes{i,3});
    
    % filter da signal for maximum determination
    sigma_filter = 2;
    data = lanes{i,2};
    threeSigma = ceil(3*sigma_filter);
    g = exp(-(-threeSigma:threeSigma).^2/2/sigma_filter^2);
    da = conv(data, g, 'same') ./ conv(ones(size(data)), g, 'same');
    
    % FIND PEAKS IN DA-channel
    peaks_da = find_peaks1d(da, 20, 0.5*max(da), 1); % window size 20, h_min=0.5*max, 1 = absolute height
 %   y_mean = peaks_da(end)+1; 
    
    y_max_max(i,:) = [lanes{i,4}(j1) lanes{i,4}(peaks_da(end)+1) lanes{i,4}(j3)];
    %y_max_max(i,:) = [lanes{i,4}(j1) lanes{i,4}(j2) lanes{i,4}(j3)];
end

subplot(2, 1, 1)
myleg = zeros(3,1);
myleg(1:3) = plot(1:size(lanes,1),y_max_max(:,1), 'g.-', 1:size(lanes,1),y_max_max(:,2), 'b.-', 1:size(lanes,1),y_max_max(:,3), 'r.-', 'MarkerSize', 15);

ylabel('Migration distance [pixel]')
set(gca, 'YDir', 'reverse')
set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',lanes(:,6),'Fontsize', 12)
set(gca, 'XLim', [0 size(lanes,1)+1])
xticklabel_rotate([1:size(lanes,1)],90,lanes(:,6))
title('Location of maxima')
legend(myleg, {'D -> D', 'D -> A', 'A -> A'}, 'location', 'best')

subplot(2, 1, 2)
plot(1:size(lanes,1),width(:,1), 'g.-', 1:size(lanes,1),width(:,2), 'b.-', 1:size(lanes,1),width(:,3), 'r.-', 'MarkerSize', 15)
xlabel('Lane')
ylabel('Width of band [pixel]')
set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',lanes(:,6),'Fontsize', 12)
set(gca, 'XLim', [0 size(lanes,1)+1])
xticklabel_rotate([1:size(lanes,1)],90,lanes(:,6))
title('Width of band from fit')


print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'MigrationDistance_Width.tif'])



%% BANDS and FRET
label_sz = 12;
close all
fig_dim =[25 20];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
E_app = [];
l = [1:size(lanes,1)];
for i=1:size(lanes,1)
    E_app = [E_app ; I_sum(i,2) ./ (I_sum(i,2) + I_sum(i,1))];
end


subplot(4, 2 , 1)
plot_subimage(dd_bg,  area,[1 5])
ylabel('D -> D')
set(gca, 'XTick', [],  'YTick', [])
hold on
for i=1:size(bands,1)
    rectangle('Position', bands{i,2}, 'EdgeColor', 'r')
end


subplot(4, 2, 2)
bar(l,I_sum(:,1), 'g')
xlabel('Lane nr.', 'Fontsize', label_sz)
%ylabel('Integrated intensity [a.u.]', 'Fontsize', label_sz)


subplot(4, 2 , 3)
plot_subimage(da_cor,  area,[1 3])
ylabel('D -> A')
set(gca, 'XTick', [],  'YTick', [])
hold on
for i=1:size(bands,1)
    rectangle('Position', bands{i,2}, 'EdgeColor', 'r')
end


subplot(4, 2, 4)
bar(l, I_sum(:,2), 'b')
xlabel('Lane nr.', 'Fontsize', label_sz)
ylabel('Integrated intensity [a.u.]', 'Fontsize', label_sz)


subplot(4, 2 , 5)
plot_subimage(aa_bg,  area,[1 3])
ylabel('A -> A')
set(gca, 'XTick', [],  'YTick', [])
hold on
for i=1:size(bands,1)
    rectangle('Position', bands{i,2}, 'EdgeColor', 'r')
end

subplot(4, 2, 6)
bar(l, I_sum(:,3), 'r')
xlabel('Lane nr.', 'Fontsize', label_sz)
%ylabel('Integrated intensity [a.u.]', 'Fontsize', label_sz)

subplot(4,2, 8)
bar(l, E_app, 'k')
xlabel('Lane nr.', 'Fontsize', label_sz)
ylabel('FRET Efficiency', 'Fontsize', label_sz)

print(cur_fig, '-dtiff', '-r500' , [path_out_plots filesep 'Integration_region.tif']); %save figure

%%
close all
fig_dim =25*[1 size(aa_bg,1)/size(aa_bg,2)];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

plot_image(aa_bg, [0.1 0.2])
%ylabel('A -> A')
set(gca, 'XTick', [],  'YTick', [])
hold on
for i=1:size(bands,1)
    rectangle('Position', bands{i,2}, 'EdgeColor', 'r')
    text(bands{i,2}(1)+bands{i,2}(3)/2, bands{i,2}(2), num2str(i), 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Center', 'Color', [1 0 0], 'FontSize', 5, 'Linewidth', 0.5)
end
print(cur_fig, '-dtiff', '-r500' , [path_out_plots filesep 'Integration_region_aa.tif']); %save figure


%% FRET FOR EACH BAND
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

E_app = [I_sum(:,2) ./ (I_sum(:,2) + I_sum(:,1))];
E = [I_sum(:,2) ./ (I_sum(:,2) + gamma*I_sum(:,1))];


subplot(2, 1, 1)
plot(1:size(lanes,1),I_sum(:,1), 'g.-', 1:size(lanes,1),I_sum(:,2), 'b.-', 1:size(lanes,1),I_sum(:,3), 'r.-', 'MarkerSize', 15)
%xlabel('Lane', 'Fontsize', label_sz)
ylabel('Integrated intensity [a.u.]', 'Fontsize', label_sz)
set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',lanes(:,6),'Fontsize', 12)
set(gca, 'XLim', [0 size(lanes,1)+1])
xticklabel_rotate([1:size(lanes,1)],90,lanes(:,6))
legend({'D->D', 'D->A', 'A->A'}, 'location', 'best')

subplot(2,1, 2)
%bar(1:size(lanes,1), E_app, 'k')
plot(1:size(lanes,1), E_app, 'k.-', 1:size(lanes,1), E, 'b.-')
%xlabel('Lane', 'Fontsize', label_sz)
ylabel('FRET Efficiency', 'Fontsize', label_sz)
set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',lanes(:,6), 'Fontsize', 12)
set(gca, 'XLim', [0 size(lanes,1)+1])
xticklabel_rotate([1:size(lanes,1)],90,lanes(:,6))
legend({'E_app, gamma = 1', ['E_true, gamma = ' num2str(round(gamma*100)/100)]}, 'location', 'best')

print(cur_fig, '-dtiff','-r600', [path_out_plots filesep 'FRET_curve.tif'])

%% FRET FOR EACH BAND
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

E_app = [I_sum(:,2) ./ (I_sum(:,2) + I_sum(:,1))];
E = [I_sum(:,2) ./ (I_sum(:,2) + gamma*I_sum(:,1))];


subplot(2, 1, 1)
plot(1:size(lanes,1),I_sum(:,1), 'g.-', 1:size(lanes,1),I_sum(:,2), 'b.-', 1:size(lanes,1),I_sum(:,3), 'r.-', 'MarkerSize', 15)
%xlabel('Lane', 'Fontsize', label_sz)
ylabel('Integrated intensity [a.u.]', 'Fontsize', label_sz)
set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',lanes(:,6),'Fontsize', 12)
set(gca, 'XLim', [0 size(lanes,1)+1])
xticklabel_rotate([1:size(lanes,1)],90,lanes(:,6))
legend({'D->D', 'D->A', 'A->A'}, 'location', 'best')

subplot(2,1, 2)
%bar(1:size(lanes,1), E_app, 'k')
plot(1:size(lanes,1), E_app, 'k.-', 1:size(lanes,1), E, 'b.-')
%xlabel('Lane', 'Fontsize', label_sz)
ylabel('FRET Efficiency', 'Fontsize', label_sz)
set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',lanes(:,6), 'Fontsize', 12)
set(gca, 'XLim', [0 size(lanes,1)+1])
set(gca, 'YLim', [0 1])
xticklabel_rotate([1:size(lanes,1)],90,lanes(:,6))

print(cur_fig, '-dtiff','-r600', [path_out_plots filesep 'FRET_curve_standard.tif'])



%% plot A->A lane for each lane
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
names = cell(1, size(lanes,1) );
cc = varycolor(size(lanes,1));
ylim = [0 0];
for i=1:size(lanes,1)
   [a b] = max(lanes{i,3});
   h(i) = plot(lanes{i,4}, lanes{i,3}, 'Color', cc(i,:));
   names{i} =  lanes{i,6};
   ylim(2) = max(ylim(2), a);
   ylim(1) = min(ylim(1), min(lanes{i,3}));
end
legend(h(1:size(lanes,1)), names{1:size(lanes,1)} )
title('A -> A Channel, raw profiles')
set(gca, 'YLim', ylim )
set(gca, 'XLim', [lanes{1,4}(1) lanes{1,4}(end)])
xlabel('Migration Distance [Pixel]')
ylabel('$I_{A \rightarrow A}$', 'Interpreter', 'Latex')
print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'AA_comparison.tif'])

%% plot A->A lane for each lane, shifted
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
names = cell(1, size(lanes,1) );
cc = varycolor(size(lanes,1));
ylim = [0 0];
for i=1:size(lanes,1)
   [a b] = max(lanes{i,3});
   h(i) = plot(lanes{i,4}-lanes{i,4}(b), lanes{i,3}, 'Color', cc(i,:));
   names{i} =  lanes{i,6};
   ylim(2) = max(ylim(2), a);
   ylim(1) = min(ylim(1), min(lanes{i,3}));
end
legend(h(1:size(lanes,1)), names{1:size(lanes,1)} )
set(gca, 'YLim', ylim )
%set(gca, 'XLim', [lanes{1,4}(1) lanes{1,4}(end)]-lanes{1,4}(b))
title('Shifted lanes')
xlabel('Migration Distance [Pixel]')
ylabel('$I_{A \rightarrow A}$', 'Interpreter', 'Latex')
print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'AA_comparison_shifted.tif'])



%% plot A->A lane for each lane, shifted and normalized
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
names = cell(1, size(lanes,1) );
cc = varycolor(size(lanes,1));
ylim = [0 0];
for i=1:size(lanes,1)
   [a b] = max(lanes{i,3});
   h(i) = plot(lanes{i,4}-lanes{i,4}(b), lanes{i,3}./a, 'Color', cc(i,:));
   names{i} =  lanes{i,6};
   ylim(2) = max(ylim(2), 1);
   ylim(1) = min(ylim(1), min(lanes{i,3}/a));
end
legend(h(1:size(lanes,1)), names{1:size(lanes,1)} )
set(gca, 'YLim', ylim )
%set(gca, 'XLim', [lanes{1,4}(1) lanes{1,4}(end)]-lanes{1,4}(b))
title('Shifted and Normalized lanes')
xlabel('Migration Distance [Pixel]')
ylabel('$I_{A \rightarrow A}$', 'Interpreter', 'Latex')
print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'AA_comparison_shifted_normalized.tif'])



%% plot A->A lane for each lane, shifted and normalized to SUM
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
names = cell(1, size(lanes,1) );
cc = varycolor(size(lanes,1));
ylim = [0 0];
for i=1:size(lanes,1)
   [a b] = max(lanes{i,3});
   s = sum(lanes{i,3});
   h(i) = plot(lanes{i,4}-lanes{i,4}(b), lanes{i,3}./s, 'Color', cc(i,:));
   names{i} =  lanes{i,6};
   ylim(2) = max(ylim(2), a./s);
   ylim(1) = min(ylim(1), min(lanes{i,3}/a));
end
legend(h(1:size(lanes,1)), names{1:size(lanes,1)} )
set(gca, 'YLim', ylim )
%set(gca, 'XLim', [lanes{1,4}(1) lanes{1,4}(end)]-lanes{1,4}(b))
title('Shifted and Normalized on sum of lane')
xlabel('Migration Distance [Pixel]')
ylabel('$I_{A \rightarrow A}$', 'Interpreter', 'Latex')
print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'AA_comparison_shifted_normalized_sum.tif'])




%% plot profiles seperately
close all
fig_dim =[20 15];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

for i=1:size(lanes,1)
    subplot(3, 1, 1)
    plot(lanes{i,4},lanes{i,1},'g'),     hold on
    vline(lanes{i,4}(sum_limits(i,1) +  sum_limits(i,2)), 'k--');
    vline(lanes{i,4}(sum_limits(i,1) -  sum_limits(i,2)), 'k--');
    legend({'D -> D'})
    ylabel('Intensity [a.u.]')
    set(gca, 'XLim', [min(lanes{i,4}) max(lanes{i,4})])
    title([lanes{i,6} ' (' num2str(i) ' of ' num2str(size(lanes, 1)) ')'])
    hold off

    subplot(3, 1, 2)
    plot( lanes{i,4},lanes{i,2},'b' ), hold on
    vline(lanes{i,4}(sum_limits(i,1) +  sum_limits(i,2)), 'k--');
    vline(lanes{i,4}(sum_limits(i,1) -  sum_limits(i,2)), 'k--');
    legend({'D -> A'})
    ylabel('Intensity [a.u.]')
    set(gca, 'XLim', [min(lanes{i,4}) max(lanes{i,4})])
    hold off

    subplot(3, 1, 3)
    plot( lanes{i,4},lanes{i,3},'r' ), hold on
    vline(lanes{i,4}(sum_limits(i,1) +  sum_limits(i,2)), 'k--');
    vline(lanes{i,4}(sum_limits(i,1) -  sum_limits(i,2)), 'k--');
    legend({'A -> A'})
    set(gca, 'XLim', [min(lanes{i,4}) max(lanes{i,4})])

    xlabel('Pixel')
    ylabel('Intensity [a.u.]')
    hold off

    print(cur_fig, '-dtiff', '-r300', [path_out_profiles filesep 'profile_seperate_lane_' sprintf('%.02i',i)])
end

close all
display('...done')

%%

E_app = [I_sum(:,2) ./ (I_sum(:,2) + I_sum(:,1))];
E = [I_sum(:,2) ./ (I_sum(:,2) + gamma*I_sum(:,1))];

l = dsDNA_distances([10:40])/10; %nm assuming B-Form DNA
E_calib = E_app(1:31);
E_gfp = E_app(32:44);
%E_gfp = E_app([32:43 45]); % gel from 2014-03-07

%
m = [
0   0   0.01
3	132	47.1241
3	157	23.3387
3	198	23.0822
3	204	27.6979
3	212	39.3466
26	132	14.4572
%26	157	35.1685
%26	198	37.1143
26	212	25.7755
132	157	39.5715
132	198	38.3533
132	204	29.3513
132	212	30.8227
157	198	15.6155
%157	204	30.9032
%198	204	20.1488
%198	212	42.9572
];
%

%{
% for gel from 2014-03-07
m = [
3	132	47.1241
0   0   0.01
3	198	23.0822
3	157	23.3387
3	212	39.3466
3	204	27.6979
26	212	25.7755
26	132	14.4572
132	198	38.3533
132	157	39.5715
132	212	30.8227
132	204	29.3513
157	198	15.6155
];
%}
d_crystal = m(:,3)/10;

gfp = cell(1,size(m,1));
for i=1:size(m,1)
    gfp{i} = [num2str(m(i,1)) '-' num2str(m(i,2))];
end

cc = varycolor(length(d_crystal));

%%
close all
fig_dim = [10 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

plot( interp1(E_calib, l, [min(E_calib):0.001:max(E_calib)], 'linear'), [min(E_calib):0.001:max(E_calib)], 'b', l, E_calib, 'k.'), hold on
legend({'Linear Interpolation', 'Calibration data'})

set(gca, 'XLim', [0 20], 'YLim', [0 1])
xlabel('C5-C5 distance [nm]'), ylabel('app. FRET Efficiency')
for i=1:length(d_crystal)
    hline(E_gfp(i), {'--', 'Color', cc(i,:), 'LineWidth', 1});
end

print(cur_fig, '-dtiff','-r500',  [path_out_plots filesep 'calibration_curve.tif'])


%%
close all
fig_dim = [20 7.5];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

d_exp =  interp1(E_calib,l, E_gfp, 'linear'); %measured distance from interpolation

cc = varycolor(length(d_crystal));

ylim = [0 1];

subplot(1, 2, 1)
plot( interp1(E_calib, l, [min(E_calib):0.001:max(E_calib)], 'linear'), [min(E_calib):0.001:max(E_calib)], 'b', l, E_calib, 'k.'), hold on
legend({'Linear Interpolation', 'Calibration data'})

set(gca, 'XLim', [0 20], 'YLim', ylim)
xlabel('C5-C5 distance [nm]'), ylabel('app. FRET Efficiency')
for i=1:length(d_crystal)
    hline(E_gfp(i), {'--', 'Color', cc(i,:), 'LineWidth', 1});
end

subplot(1, 2, 2)
h = zeros(length(d_crystal), 1);
for i=1:length(d_crystal)
    bar(i, E_gfp(i), 'FaceColor', cc(i,:)), hold on
end
for i=1:length(d_crystal)
    h(i) =    hline(E_gfp(i), {'--', 'Color', cc(i,:), 'LineWidth', 1});
end

ylabel('app. FRET Efficiency')
xlabel('GFP-mutant')
set(gca, 'YLim', ylim)

print(cur_fig, '-dtiff','-r500',  [path_out_plots filesep 'calibration_samples.tif'])

%%
close all
fig_dim = [15 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
h = zeros(length(d_exp)+1,1)' ;
for i=[1:length(d_exp)]
    h(i) = plot(d_crystal(i), d_exp(i), 'o', 'MarkerFaceColor', cc(i,:), 'MarkerEdgeColor', [0 0 0], 'Markersize', 10); hold on
 end
set(gca, 'XLim', [0 14], 'YLim', [0 14]), axis square
set(gca, 'XTick', [0:2:14], 'YTick', [0:2:14])

i=[1:length(d_exp)];
[c, S] = polyfit(d_crystal(i), d_exp(i), 1);

ste = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df);


h(end) = plot([0 15], c(1)*[0 15]+c(2));

legend(h, [gfp {['Fit, slope = ' num2str(round(c(1)*100)/100) '+-' num2str(round(ste(1)*100)/100) ', offset = ' num2str(round(c(2)*100)/100)  '+-' num2str(round(ste(2)*100)/100)  ]}], 'Location', 'Southeast')
grid on
xlabel('Separation from crystal structure [nm]')
ylabel('Separation measured [nm]')

print(cur_fig, '-dtiff','-r500',  [path_out_plots filesep 'experimental_vs_crystal.tif'])

%%
close all
fig_dim = [20 20];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
h = zeros(length(d_exp)+1,1)' ;
for i=[1:length(d_exp)]
    h(i) = plot(d_crystal(i), d_exp(i)-c(2), 'o', 'MarkerFaceColor', cc(i,:), 'MarkerEdgeColor', [0 0 0], 'Markersize', 10); hold on
end
set(gca, 'XLim', [0 6], 'YLim', [0 6]), axis square
set(gca, 'XTick', [0:2:14], 'YTick', [0:2:14])
h(end) = plot([0 15], c(1)*[0 15]);

legend(h, [gfp {['Fit, slope = ' num2str(round(c(1)*100)/100) '+-' num2str(round(ste(1)*100)/100) ', offset = ' num2str(round(c(2)*100)/100)  '+-' num2str(round(ste(2)*100)/100)  ]}], 'Location', 'best')
grid on
xlabel('Separation from crystal structure [nm]')
ylabel('Separation measured [nm]')

print(cur_fig, '-dtiff','-r500',  [path_out_plots filesep 'experimental_vs_crystal_offset-corrected.tif'])

%%
dlmwrite([path_out_plots filesep 'mutant_distance_data.txt'] , [m d_exp*0.34], '\t')
dlmwrite([path_out_plots filesep 'calibration_data.txt'] , [l E_calib], '\t')

