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
correction = 1;
leak_dir  = calculate_corrections(dd_bg, da_bg, aa_bg, [path_out filesep prefix_out '_correction.txt']);
da_cor = da_bg - leak_dir(1,1).*dd_bg - leak_dir(2,1).*aa_bg;%-leak_dir(1,2)-leak_dir(2,2); 


%% gamma correction
gamma = 1;

%% find lanes
[auto_pos , area] = find_lanes(aa_bg+dd_bg);
area = [min(area(:,1)) min(area(:,2)) max(area(:,1)+area(:,3))-min(area(:,1)) max(area(:,2)+area(:,4))-min(area(:,2))]; %area which surrounds all subareas
n_lanes = size(auto_pos,1);

%% find control lanes
button = questdlg('Select controls?','Controls','Yes','No', 'Yes');
correction = 1;
if strcmp(button, 'Yes') 
    [auto_pos_control , area_control] = find_lanes(aa_bg+dd_bg);
    auto_pos = [auto_pos; auto_pos_control]; % append to auto_pos
    n_lanes = n_lanes + size(auto_pos_control,1);
end

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


%% fit gaussian close to maxima
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
   

    % FIND PEAKS IN DA-channel
%    peaks_da = find_peaks1d(da, 20, 0.5*max(da), 1); % window size 20, h_min=0.5*max, 1 = absolute height
%    y_mean = peaks_da(end);  % use the last found peak (leading band)
%    dy = 5; % integration width
    
    % get the integration region from maximum
    [dd_max, dd_imax] = max(lanes{i,1});
    [da_max, da_imax] = max(lanes{i,2});
    [aa_max, aa_imax] = max(lanes{i,3});
    
    
    y_weighted = round(     (dd_max*dd_imax + aa_max*aa_imax)/(dd_max+aa_max)    );
    
    y_mean = y_weighted; %round((dd_imax + aa_imax)/2);
    dy = 5;
    
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
    
    %{
   y = lanes{i,4};
   plot(lanes{i,4},lanes{i,1},'g',  lanes{i,4},lanes{i,2},'b',  lanes{i,4},lanes{i,3},'r' ), hold on
   plot(y, gauss1d(p_dd, y), 'g--', y, gauss1d(p_da, y), 'b--', y, gauss1d(p_aa, y), 'r--')
   vline(y(y_mean-dy), 'k--')
   vline(y(y_mean+dy), 'k--')
   vline(y(y_mean), 'k')
   vline(y(y_weighted), 'b')

   hold off
   pause
   %pause(1)
    %}
end
close all

%% calculate ratios
ratio = zeros(n_lanes, 3);

for i=1:n_lanes
    pos = bands{i,2};
    DD = dd_bg( pos(2):pos(2)+pos(4) , pos(1):pos(1)+pos(3) );
    AA = aa_bg( pos(2):pos(2)+pos(4) , pos(1):pos(1)+pos(3) );
    DA = da_cor( pos(2):pos(2)+pos(4) , pos(1):pos(1)+pos(3) );
    
    ratio(i,1) = calculateRation(DD, AA, 0); % DD / AA
    ratio(i,2) = calculateRation(AA, DD, 0); % AA / DD
    ratio(i,3) = calculateRation(DD, DA, 0); % DD / DA
end

%% save all the data
disp('Saving data...')
save([path_out filesep prefix_out '_data.mat'])

%% sava data for fret only
save([path_out filesep prefix_out '_data_ratio.mat'], 'ratio','I_sum')

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
set(gca, 'YLim', [0 20])
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




%% plot D->D lane for each lane
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
names = cell(1, size(lanes,1) );
cc = varycolor(size(lanes,1));
ylim = [0 0];
for i=1:size(lanes,1)
   [a b] = max(lanes{i,1});
   h(i) = plot(lanes{i,4}, lanes{i,1}, 'Color', cc(i,:));
   names{i} =  lanes{i,6};
   ylim(2) = max(ylim(2), a);
   ylim(1) = min(ylim(1), min(lanes{i,1}));
end
legend(h(1:size(lanes,1)), names{1:size(lanes,1)} )
title('D -> D Channel, raw profiles')
set(gca, 'YLim', ylim )
set(gca, 'XLim', [lanes{1,4}(1) lanes{1,4}(end)])
xlabel('Migration Distance [Pixel]')
ylabel('D -> D')
print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'DD_comparison.tif'])

%% plot D->D lane, shifted
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
names = cell(1, size(lanes,1) );
cc = varycolor(size(lanes,1));
ylim = [0 0];
for i=1:size(lanes,1)
   [a b] = max(lanes{i,1});
   h(i) = plot(lanes{i,4}-lanes{i,4}(b), lanes{i,1}, 'Color', cc(i,:));
   names{i} =  lanes{i,6};
   ylim(2) = max(ylim(2), a);
   ylim(1) = min(ylim(1), min(lanes{i,1}));
end
legend(h(1:size(lanes,1)), names{1:size(lanes,1)})
title('D -> D Channel, raw profiles')
set(gca, 'YLim', ylim )
%set(gca, 'XLim', [lanes{1,4}(1) lanes{1,4}(end)])
xlabel('Migration Distance [Pixel]')
ylabel('D -> D')
print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'DD_comparison_shifted.tif'])

%% plot D->D lane, shifted and normamlized
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
h = zeros(1, size(lanes,1));
names = cell(1, size(lanes,1) );
cc = varycolor(size(lanes,1));
ylim = [0 0];
for i=1:size(lanes,1)
   [a b] = max(lanes{i,1});
   h(i) = plot(lanes{i,4}-lanes{i,4}(b), lanes{i,1}./a, 'Color', cc(i,:));
   names{i} =  lanes{i,6};
   ylim(2) = max(ylim(2), 1);
   ylim(1) = min(ylim(1), min(lanes{i,1}./a));
end
legend(h(1:size(lanes,1)), names{1:size(lanes,1)})
title('D -> D Channel, raw profiles')
set(gca, 'YLim', ylim )
%set(gca, 'XLim', [lanes{1,4}(1) lanes{1,4}(end)])
xlabel('Migration Distance [Pixel]')
ylabel('D -> D')
print(cur_fig, '-dtiff', '-r500', [path_out_plots filesep 'DD_comparison_shifted_normalized.tif'])





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

