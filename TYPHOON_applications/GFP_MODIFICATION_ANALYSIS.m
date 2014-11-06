%%
clear all, close all, clc
run('my_prefs')
path0 = cd;

%% load image
%nimg_input = inputdlg({'Number of images:'}, 'NUmber of images' , 1, {'1'} );
n_img = 2;%str2double(nimg_input{1});

filenames = cell(n_img, 1);
pathnames = cell(n_img, 1);

last_dir = data_dir;
for i=1:n_img
    cd(last_dir)
    if i==1
        [filenames{i} pathnames{i}]=uigetfile('*.tif','Select [Cy2] image:');
    end
    if i==2
        [filenames{i} pathnames{i}]=uigetfile('*.tif','Select [Cy5] image:');
    end
    last_dir = pathnames{i};
end
cd(path0)

%% create output folder
pname = inputdlg({'Output folder and prefix:'}, 'Output folder and prefix' , 1, {[filenames{1}(1:size(filenames{1},2)-10) '_analysis']} );
prefix_out = pname{1};
path_out = [pathnames{1} prefix_out ];
mkdir(path_out);

%% load and bg correct images
images = cell(n_img, 1);
img_bg = cell(n_img, 1);
for i=1:n_img
    images{i} = double(imread([pathnames{i} filesep filenames{i}]));  %load
    plot_image_ui(images{i})
    button = questdlg('Rotate?','Rotate','Rotate','No','No');
    if strcmp(button,'Rotate') %load old data
        images{i} = imrotate(images{i}, -90);
    end
    close all
    img_bg{i} = bg_correct_ui(images{i}, 'Background correction');    %bg correct
    close all    
end


%% select bands by hand
close all
subplot(1, 2, 1)
plot_image(img_bg{1}, 0.1)
subplot(1, 2, 2)
plot_image(img_bg{2}, 0.1)
options.WindowStyle='normal';
prompt={'How many bands'};
def={'4'};
tmp = inputdlg(prompt, 'How many bands', 1, def, options);
n_bands = str2double(tmp(1));
close all
%%
I = zeros(n_bands, n_img);
I_max = zeros(n_bands, n_img);
areas = zeros(n_bands, 4);

[I, areas] = integrate_areas(img_bg, n_bands, 1);

for i=1:n_bands
    pos = areas(i,:) ;
    for j = 1:n_img
        I_max(i,j) = max(max(img_bg{j}( pos(2):pos(2)+pos(4)  , pos(1):pos(1)+pos(3)  )));
    end
end

%% assining names
names = cell(n_bands,1);
names2 = cell(n_bands/2,1);
for i=1:2:n_bands
    band_name = inputdlg({'Name of sample:'}, 'Sample names' , 1, {['Sample ' num2str((i+1)/2)]} );
    names{i} = [band_name{1} ' 1xlabeled'];
    names{i+1} = [band_name{1} ' 2xlabeled'];
    names2{(i+1)/2} = band_name{1} ;


    
end




%%
%colormap('gray'), colorbar

%imagesc(img_bg)
for i=1:n_img
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual');
    colormap('gray'), colorbar

    plot_image(img_bg{i}, 0.1)
    title(filenames{i})
    for j=1:n_bands
       rectangle('Position', areas(j,:), 'EdgeColor', 'r')
       text(double(areas(j,1)), double(areas(j,2)), num2str(j), 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'top')
    end 
    %plot_image_ui(img_bg{i})

    
    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.eps']); %save figure
    print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.png']); %save figure
    
    close all
end

%% find lanes

[auto_pos , area] = find_lanes(sum(sum(img_bg{1})).*(img_bg{1}./sum(sum(img_bg{1}))+img_bg{2}./sum(sum(img_bg{2}))));
%%
n_lanes = size(auto_pos,1);
lanes = cell(n_lanes, n_img+1);

for i=1:n_lanes
    lanes{i,1} = (auto_pos(i,2):auto_pos(i,2)+auto_pos(i,4))'; %y
    for j=1:n_img
        cur_lane = img_bg{j}(auto_pos(i,2):auto_pos(i,2)+auto_pos(i,4)  ,   auto_pos(i,1):auto_pos(i,1)+auto_pos(i,3));
        lanes{i,1+j} = transpose( sum(transpose( cur_lane) ) ); % profiles
    end
end

%%
save([path_out filesep prefix_out '_data'] )



%%
%colormap('gray'), colorbar

%imagesc(img_bg)
for i=1:n_img
 %   cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual');
 %   colormap('gray'), colorbar

    %plot_image(img_bg{i},filenames{i}, 12 )
    plot_image_ui(img_bg{i})
    for j=1:n_bands
       rectangle('Position', areas(j,:), 'EdgeColor', 'r')
       text(double(areas(j,1)), double(areas(j,2)), num2str(j), 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'top')
    end 
    %plot_image_ui(img_bg{i})

    
   % print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.eps']); %save figure
   % print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.png']); %save figure
    pause
    close all
end


%% PLOTTING STUFF
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I)
set(gca, 'XTick',1:n_bands, 'XTickLabel',names,'Fontsize', 12)
xticklabel_rotate([1:n_bands],90,names)

set(gca, 'XTick', 1:n_bands)
ylabel('Intensity')
legend({'GFP', 'DNA'})
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_IntesityRaw.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_IntesityRaw.png']); %save figure

%% PLOTTING STUFF
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I(1:2:n_bands,:)+I(2:2:n_bands,:))

set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands)
ylabel('Sum of labeled- and unlabeled-band')
%ylabel('Sum of monomer and dimer band')
legend({'GFP', 'DNA'})

print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_sum.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_sum.png']); %save figure

%%
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


A = [I(2:2:n_bands,1)./(I(1:2:n_bands,1)+I(2:2:n_bands,1)) I(2:2:n_bands,2)./(2*I(1:2:n_bands,2)+I(2:2:n_bands,2))];
bar(A)
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands/2)
%set(gca, 'YLim', [0 1])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
%ylabel('[dimer] / ([dimer] + [monomer])')
ylabel('[GFP-2DNA] / ([GFP-2DNA] + [GFP-1DNA])')

legend({'GFP-channel', 'DNA-channel'})

print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_normalizedIntensities.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_normalizedIntensities.png']); %save figure


%%
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

yield = I(2:2:n_bands,1)./(I(1:2:n_bands,1)+I(2:2:n_bands,1));
bar(yield, 'g');
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands/2)
set(gca, 'YLim', [0 1.2])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('[GFP+2DNA] / ([GFP+1DNA] + [GFP+2DNA])')

legend({'GFP-channel'})

for i=1:size(yield,1)
    text(i,0.5, [num2str(round(yield(i)*100)) ' %'], 'Color', 'black', 'FontSize', 15, 'HorizontalAlignment','center')
end

print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_normalizedIntensities2.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_normalizedIntensities2.png']); %save figure


for i=1:size(yield,1)
    display(['Yield of ' names2{i} ' from GFP-channel: ' num2str( round(yield(i)*100)  ) ' %'])
end


%%
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

yield2 = I(2:2:n_bands,2)./(I(1:2:n_bands,2)*2+I(2:2:n_bands,2));
bar(yield2, 'r');
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands/2)
set(gca, 'YLim', [0 1.2])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('[GFP+2DNA] / ([GFP+1DNA] + [GFP+2DNA])')

legend({'Reporter-channel'})

for i=1:size(yield2,1)
    text(i,0.5, [num2str(round(yield2(i)*100)) ' %'], 'Color', 'black', 'FontSize', 15, 'HorizontalAlignment','center')
end

print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_normalizedIntensities3.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_normalizedIntensities3.png']); %save figure


for i=1:size(yield,1)
    display(['Yield of ' names2{i} ' from Reporter-channel: ' num2str( round(yield2(i)*100)  ) ' %'])
end


%% plot reporter to gfp ratio
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(1:n_bands, I(1:1:n_bands,2) ./ I(1:n_bands,1));
set(gca, 'XTick',1:n_bands, 'XTickLabel',names(1:end),'Fontsize', 12)
xticklabel_rotate([1:n_bands],90,names(1:end))

set(gca, 'XTick', 1:n_bands)
%set(gca, 'YLim', [0 2])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('[Reporter] / [GFP]')
title('Ratio of Reporter- to GFP-Signal')

for i=1:n_bands
    text(i,0.5, [num2str(round(10*I(i,2) ./ I(i,1))./10)], 'Color', 'red', 'FontSize', 14, 'HorizontalAlignment','center')
end


print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_relIntensity.png']); %save figure



%% write rgb tif

rgb = write2rgb(img_bg{2}, img_bg{1}, zeros(size(img_bg{1})), [path_out filesep prefix_out '_rgb.tif'] );
rgb_scaled = write2rgb((2^16-1)*img_bg{2}/max(I_max(:,2)), (2^16-1)*img_bg{1}/max(I_max(:,1)), zeros(size(img_bg{1})), [path_out filesep prefix_out '_rgb_scaled.tif'] );


%%
colors = {'g', 'r'};
channel = {'GFP', 'Reporter'};

close all

fig_dim =[13 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


colormap gray
for i=1:n_lanes
    for j=1:n_img
        subplot(1, n_img*2, (j-1)*2+1)
        imagesc(img_bg{j}(auto_pos(i,2):auto_pos(i,2)+auto_pos(i,4) , auto_pos(i,1):auto_pos(i,1)+auto_pos(i,3)))
       set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        title(names2{i})
        
        subplot(1, n_img*2, (j-1)*2+2)    
        plot(lanes{i,1+j}, lanes{i,1},colors{j} )
        set(gca,'YDir','reverse');
        set(gca, 'YLim', [lanes{i,1}(1) lanes{i,1}(end) ])
        
        scale = [min(lanes{i,1+j}) max(lanes{i,1+j}) ];
        ds = scale(2)-scale(1);
        scale(1) = scale(1)-0.1*ds; scale(2) = scale(2)+0.1*ds; 

        
        set(gca, 'XLim', scale)
        legend(channel{j})
    end

    
    print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_lane_' num2str(i) '.png']); %save figure

    
end


%% plot all lanes in one image
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


subplot(2,1,1)
h = zeros(n_lanes,1);
hold all
for i=1:n_lanes
   h(i) = plot(lanes{i,1}, lanes{i,2});
end
legend(h, names2)
set(gca, 'XLim', [lanes{1,1}(1) lanes{1,1}(end)])
title('GFP channel')

subplot(2,1,2)
h = zeros(n_lanes,1);
hold all
for i=1:n_lanes
   h(i) = plot(lanes{i,1}, lanes{i,3});
end
legend(h, names2)
set(gca, 'XLim', [lanes{1,1}(1) lanes{1,1}(end)])
title('Reporter channel')


print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_all_lanes.png']); %save figure


%% plot relative intensity
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

for i=1:n_lanes
    
    subplot(2, 1, 1)
    plot(lanes{i,1}, lanes{i,2}, 'g-', lanes{i,1}, lanes{i,3}, 'r-')
    set(gca, 'XLim', [lanes{i,1}(1) lanes{i,1}(end)])
    scale = [min([lanes{i,2}; lanes{i,3}]) max([lanes{i,2}; lanes{i,3}])  ];
    ds = scale(2)-scale(1);
    scale(1) = scale(1)-0.1*ds; scale(2) = scale(2)+0.1*ds; 
    set(gca, 'YLim', scale)
    legend({'GFP', 'Reporter'})
    title(names2{i})
    
    subplot(2, 1, 2)
    plot(lanes{i,1}, lanes{i,3}./lanes{i,2}, 'b-')
    set(gca, 'XLim', [lanes{i,1}(1) lanes{i,1}(end)])
    set(gca, 'YLim', [0 4])
    legend({'Reporter / GFP'})
    
    print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_lane_relIntensity_' num2str(i) '.png']); %save figure

    
end

%%
display('DONE')

%%
button = questdlg('Use a gel-standard?','Gel standard','Yes','No','No');
if strcmp(button,'Yes') 
    tmp = inputdlg({'How many references:'}, 'How many references', 1, {'5'}, options);
    n_ref = str2double(tmp(1));
    
    [I_ref, areas_ref] = integrate_areas(img_bg(2), n_ref, 1);
    
    input = cell( 1, n_ref);
    def = cell(1, n_ref);
    for i=1:n_ref
       input{i} = ['Initial conc. of band ' num2str(i) ' [nM]'];
       def{i} = num2str((i-1)*50);
    end
    c_ref = str2double(inputdlg(input, 'Enter conc', 1, def));
end

%% fit gel standard
close all
fig_dim =[15 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
c_plot = 0:max(c_ref)/length(c_ref)/10:max(c_ref)*1.2;
cf = polyfit(c_ref, I_ref, 1);
plot(c_ref, I_ref, 'b.', 'Markersize', 15), hold on
plot(c_plot, cf(2)+cf(1)*c_plot, 'r--')
hline(I(2:2:end,2), 'b-')
xlabel('Concentration [nM]'), ylabel('Intensity [a.u.]')

print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_calibration.tif']); %save figure
dlmwrite([path_out filesep prefix_out '_concentration-standard_coefficients.txt'], cf,'delimiter' , '\t')

%% estimate concentration of double-labeled species
c_final = zeros(size(I,1),1);

%display('Final concentration of double-labeled-GFP [uM]:')
c_final(2:2:end) = (I(2:2:end,2)-cf(2))/cf(1)/2; % divide by 2, since 2xlabeled species has twice the signal per GFP-molecule
c_final(1:2:end) = (I(1:2:end,2)-cf(2))/cf(1);


display('Initial concentration of double-labeled-GFP [nM]:')
c_initial = c_final*10/3

display('Initial concentration of double-labeled-GFP [nM]:')
round(c_final*10/3)
%%
%write output file
fileID = fopen([path_out filesep prefix_out '_concentrations.txt'],'w');
fprintf(fileID,'#mutant\t c_final [nM]\t c_init [nM]\n');

for i=1:size(c_final,1)
    fprintf(fileID,'%s\t%i\t%i\n',names{i}, round(c_final(i)), round(c_initial(i)));
end
fclose(fileID);
%
save([path_out filesep prefix_out '_data'] )

%%

display('Final concentration of single-labeled-GFP [uM]:')
c_final_2 = (I(1:1:end,2)-cf(2))/cf(1)/2

display('Initial concentration of single-labeled-GFP [nM]:')
c_initial = round(c_final_2*10/3*1000)
