%%
clear all
close all
clc
run('my_prefs')
path0 = cd;

%% load images
options.WindowStyle='normal';
prompt={'How many images do you want to load:'};
def={'2'};
tmp = inputdlg(prompt, 'How many images do you want to load', 1, def, options);
n_img = str2double(tmp(1));

filename = cell(n_img, 1);
pathname = cell(n_img, 1);
cd(data_dir)
for i=1:n_img
    [filename{i}, pathname{i}]=uigetfile('*.tif',['Select image ' num2str(i) ' of ' num2str(n_img)]);
    cd(pathname{i})
end
cd(path0)


%% create output folder
pname = inputdlg({'Output folder and prefix:'}, 'Output folder and prefix' , 1, {filename{1}(1:size(filename{1},2)-4)} );
prefix_out = pname{1};
path_out = [pathname{1} prefix_out ];
mkdir(path_out);
%% load, bg correct images and find lanes
img = cell(n_img, 1);
img_bg = cell(n_img, 1);

for i=1:n_img
    %bg correct
    img{i} = double(imread([pathname{i} filesep filename{i}])); 
    
    plot_image_ui(img{i})
    button = questdlg('Rotate?','Rotate','Rotate','No','No');
    if strcmp(button,'Rotate') %load old data
        img{i} = imrotate(img{i}, -90);
    end
    close all
    
    
    img_bg{i} = bg_correct_ui(img{i}, 'Background correction');
    
   
end
 %% find lanes
img_sum = img_bg{1};% ./ (sum(sum(img_bg{1}))); % max(max(img_bg{1})); %

for i=2:n_img
  img_sum = img_sum + img_bg{i};% ./ (sum(sum(img_bg{i})));  %max(max(img_bg{1}));% 
end
 
[auto_pos , area] = find_lanes(img_sum);

close all



%%
close all
n_lanes = size(auto_pos,1);


lanes = cell(n_lanes, n_img+1);

for i=1:n_lanes

    lanes{i,2} = (auto_pos(i,2):auto_pos(i,2)+auto_pos(i,4))'; %y
    for j=1:n_img
        cur_lane = img_bg{j}(auto_pos(i,2):auto_pos(i,2)+auto_pos(i,4)  ,   auto_pos(i,1):auto_pos(i,1)+auto_pos(i,3));
        lanes{i,2+j} = transpose( sum(transpose( cur_lane) ) ); % profiles
    end

    if i>1
        subplot(2, 1, 2)
        plot_image(img_bg{j}, 0.1)
        hold on
        rectangle('Position', auto_pos(i,:), 'EdgeColor', 'r')
        text(double(auto_pos(i,1)), double(auto_pos(i,2)), lanes{i-1,1}, 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'bottom')
    end

    subplot(2, 1, 1)
    plot_image(img_bg{1}, 0.1)
    hold on
    rectangle('Position', auto_pos(i,:), 'EdgeColor', 'r')
    hold off

    lane_name = inputdlg({'Name of lane:'}, 'Lane properties' , 1, {['Lane ' num2str(i)]} );
    lanes{i,1} = lane_name{1};
end
close all

%%
save([path_out filesep prefix_out '_data'] )
display('done')

%%

A_out = zeros(size(lanes{1,2},1), 1+size(lanes,1));
A_out(:,1) = lanes{1,2};

for i=1:size(lanes,1)
    A_out(:,1+i)= lanes{i,3};
end
dlmwrite([path_out filesep 'data.txt'], A_out, 'delimiter' , '\t')

%%

A_out = zeros(size(lanes{1,2},1), 1+size(lanes,1));
A_out(:,1) = lanes{1,2};

for i=1:size(lanes,1)
    A_out(:,1+i)= lanes{i,3}./auto_pos(i,3);
end
dlmwrite([path_out filesep 'data_normalized.txt'], A_out, 'delimiter' , '\t')
dlmwrite([path_out filesep 'width.txt'], auto_pos(:,3), 'delimiter' , '\t')


%% plot lanes

for i=1:n_lanes
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

    
    title(lanes{i,1})
    plot(lanes{i,2}, lanes{i,3},'b' )
    legend(lanes{i,1})
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
    set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])


    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_Intesity_' num2str(i) '.eps']); %save figure
    print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_Intesity_' num2str(i) '.png']); %save figure

    close all
    
end
    

%% plot all profiles of each image

% raw data
for j=1:n_img
    close all
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    hold all
    h = zeros(1, n_lanes);
    cc = varycolor(n_lanes);

    for i=1:n_lanes
        h(i)=plot(lanes{i,2}, lanes{i,2+j} , 'Color', cc(i,:));
    end

    legend(h,lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
    set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])


    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity'  '.png']); %save figure
end

% normalized data
for j=1:n_img
    close all
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    hold all
    h = zeros(1, n_lanes);
    cc = varycolor(n_lanes);

    for i=1:n_lanes
        h(i)=plot(lanes{i,2}, lanes{i,2+j}./sum(lanes{i,2+j}) , 'Color', cc(i,:));
    end

    legend(h,lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
    set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])


    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_normalized'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_normalized'  '.png']); %save figure
end

% normalized and shifted data
for j=1:n_img
    close all
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    hold all
    h = zeros(1, n_lanes);
    cc = varycolor(n_lanes);

    for i=1:n_lanes
        [I_max i_max] = max(lanes{i,2+j});
        p =find_peaks1d(lanes{i,2+j}, 10, max(lanes{i,2+j})/4);

        i_max = p(1);
        h(i)=plot(lanes{i,2}- lanes{i,2}(i_max), lanes{i,2+j}./sum(lanes{i,2+j}), 'Color', cc(i,:) );
    end

    legend(h,lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
    %set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])


    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_normalized_shifted'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_normalized_shifted'  '.png']); %save figure
end

%% all in one
for j=1:n_img
    close all
    fig_dim =[20 20];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    h = zeros(3, n_lanes);
    %cc = varycolor(n_lanes);
    cc = [1 0 0; 0 0 1]
    for i=1:n_lanes
        [I_max i_max] = max(lanes{i,2+j});
        p =find_peaks1d(lanes{i,2+j}, 10, max(lanes{i,2+j})/4);

        i_max = p(1);
        
        subplot(3, 1, 1)
        title('Raw data')
        h(1,i)=plot(lanes{i,2}, lanes{i,2+j}, 'Color', cc(i,:) ); hold on;
        
        subplot(3, 1, 2)
        title('Normalized to sum of whole lane')
        h(2,i)=plot(lanes{i,2}, lanes{i,2+j}./sum(lanes{i,2+j}), 'Color', cc(i,:) ); hold on;
        
        subplot(3, 1, 3)
        title('Normalized and shifted to maxima')
        h(3,i)=plot(lanes{i,2}- lanes{i,2}(i_max), lanes{i,2+j}./sum(lanes{i,2+j}), 'Color', cc(i,:) ); hold on;
    end
    subplot(3, 1, 1)
    legend(h(1,:),lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]'), ylabel('Intensity')
    subplot(3, 1, 1)
    legend(h(2,:),lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]'), ylabel('Intensity')
    subplot(3, 1, 1)
    legend(h(3,:),lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]'), ylabel('Intensity')

    %set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])


    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_all'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_all'  '.png']); %save figure
end




%%

for i=1:n_lanes
    close all
    hold all
    for j=1:n_img
        plot(lanes{i,2}, lanes{i,2+j})
    end
    pause
    
    
    
end



%% fit two gaussians

    
    closed = zeros(n_lanes,3);
    open = zeros(n_lanes,3);
%%
    
for j=1:n_img
    close all
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    hold all
    h = zeros(1, n_lanes);
    cc = varycolor(n_lanes);
    for i=1:10%n_lanes
        if i == 1 || i == 11% just one peak
           h(i)=plot(lanes{i,2}, lanes{i,2+j} , 'Color', cc(i,:));
           [p p_err] = fit_peak(lanes{i,2}, lanes{i,2+j});   
           open(i,:) = p;
           closed(i,:) = [p(1) 0 0];
        else
            h(i)=plot(lanes{i,2}, lanes{i,2+j} , 'Color', cc(i,:));
            [pfast pslow pfast_err pslow_err] = fit_2peaks(lanes{i,2}, lanes{i,2+j});
            open(i,:) = pslow;
            closed(i,:) = pfast;
        end
        pause
        
    end
end
%%
dlmwrite([path_out filesep prefix_out '_open.txt'], open, 'delimiter', '\t')
dlmwrite([path_out filesep prefix_out '_closed.txt'], closed, 'delimiter', '\t')

%%
x = closed(:,2).*closed(:,3) ./ (closed(:,2).*closed(:,3) + open(:,2).*open(:,3)) ;
ratio = [0 0.4:0.1:1 1 2]';
plot(ratio, x(1:10), 'b.-', 'markersize', 15), hold on
plot(ratio, x(11:20), 'r.-', 'markersize', 15), hold on
dlmwrite([path_out filesep prefix_out '_fraction_closed.txt'], [[ratio; ratio] x], 'delimiter', '\t')

%%
close all
  fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
   
semilogx(ratio, x(1:10), 'b.-', 'markersize', 15), hold on
semilogx(ratio, x(11:20), 'r.-', 'markersize', 15)
set(gca, 'XLim', [0.1 3])
legend({'ssLoop', 'v2Loop'}, 'location', 'northwest')
xlabel('Ratio of [Linker]/[FS]')
ylabel({'Fraction of closed structures', '[closed] / ([closed]+[open])'})
set(gca, 'YLim', [0 1])

    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_fraction'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_fraction'  '.png']); %save figure


%%

for j=1:n_img
    close all
    fig_dim =[30 7.5];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    h = zeros(1, n_lanes);
    cc = varycolor(n_lanes/2);
    cc = [cc ; cc];
    
    subplot(1, 2, 1)
    hold all
    for i=1:n_lanes/2 
        h(i)=plot(lanes{i,2}, lanes{i,2+j}./sum(lanes{i,2+j}) , 'Color', cc(i,:));
    end
    legend(h(1:10),lanes{1:10,1}, 'location', 'northeast')
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
    set(gca, 'XLim', [lanes{1,2}(1) lanes{1,2}(end)])
    set(gca, 'YLim', [0 0.05])

    
    subplot(1, 2, 2)  
    hold all
    for i=11:n_lanes
        h(i)=plot(lanes{i,2}, lanes{i,2+j}./sum(lanes{i,2+j}) , 'Color', cc(i,:));
    end
    legend(h(11:20),lanes{11:20,1}, 'location', 'northeast')
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
    set(gca, 'XLim', [lanes{1,2}(1) lanes{1,2}(end)])
    set(gca, 'YLim', [0 0.05])

    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_normalized'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_normalized'  '.png']); %save figure
end

    %%
    
    
for j=1:n_img

    close all
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    hold all
    h = zeros(1, n_lanes);
    cc = varycolor(n_lanes);

    for i=1:n_lanes
        [I_max i_max] = max(lanes{i,2+j});
        p =find_peaks1d(lanes{i,2+j}, 10, max(lanes{i,2+j})/4);



        i_max = p(1);
        h(i)=plot(lanes{i,2}- lanes{i,2}(i_max), lanes{i,2+j}./sum(lanes{i,2+j}), 'Color', cc(i,:) );
    end

    legend(h,lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]')
    ylabel('Normalized Intensity [a.u.]')
    %set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])


    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_shifted'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_' num2str(j)  '_Intesity_shifted'  '.png']); %save figure

end    
%%
for j=1:n_img

    close all
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    hold all
    h = zeros(1, n_lanes);
    cc = varycolor(n_lanes);
    for i=1:n_lanes
        [I_max i_max] = max(lanes{i,2+j});

        h(i)=plot(lanes{i,2}-lanes{i,2}(i_max), lanes{i,2+j}./I_max, 'Color', cc(i,:) );


    end

    legend(h,lanes{:,1}, 'location', 'best')
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
    %set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])


    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_img_' num2str(j) '_Intesity_normalized'  '.eps']); %save figure
    print(cur_fig, '-dpng','-r600' , [path_out filesep prefix_out '_img_' num2str(j)  '_Intesity_normalized'  '.png']); %save figure

end    
%%



for i=1:n_lanes
    fig_dim =[20 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

    
    subplot(2, 1, 1)
    plot(lanes{i,2}, lanes{i,3},'g' , lanes{i,2}, lanes{i,4},'r' )
    legend({ 'GFP', 'Reporter'})
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
                set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])
    title(lanes{i,1})

    subplot(2, 1, 2)
    plot(lanes{i,2}, lanes{i,4}./(lanes{i,3}) )
    legend({'Reporter / GFP '})
    xlabel('Migration distance [pixel]')
    ylabel('rel. Intensity')
            set(gca, 'YLim', [0 10])
            set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])

    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_relIntesity_' num2str(i) '.eps']); %save figure
    
    
    close all
    
end



%%


for i=1:n_lanes
    fig_dim =3*[10 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

    
    subplot(2, 1, 1)
    plot(lanes{i,2}, lanes{i,3},'r' , lanes{i,2}, lanes{i,4},'g' )
    legend({ 'GFP', 'Reporter'})
    xlabel('Migration distance [pixel]')
    ylabel('Intensity')
                set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])
    title(lanes{i,1})

    subplot(2, 1, 2)
    plot(lanes{i,2}, (lanes{i,4})./(lanes{i,3}) )
    legend({'GFP / FS '})
    xlabel('Migration distance [pixel]')
    ylabel('rel. Intensity')
            set(gca, 'YLim', [0 2])
            set(gca, 'XLim', [lanes{i,2}(1) lanes{i,2}(end)])

    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_relIntesity_' num2str(i) '.eps']); %save figure
    pause
    
    close all
    
end

%%
colors = {'g', 'r'};
channel = {'GFP', 'Reporter'};

close all

fig_dim =[20 20];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


colormap gray
for i=1:n_lanes
    
    for j=1:n_img
        subplot(1, n_img*2, (j-1)*2+1)
        imagesc(img_bg{j}(auto_pos(i,2):auto_pos(i,2)+auto_pos(i,4) , auto_pos(i,1):auto_pos(i,1)+auto_pos(i,3)))
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        title(lanes{i,1})
        
        subplot(1, n_img*2, (j-1)*2+2)    
        plot(lanes{i,2+j}, lanes{i,2},colors{j} )
        set(gca,'YDir','reverse');
        set(gca, 'YLim', [lanes{i,2}(1) lanes{i,2}(end) ])
        
        scale = [min(lanes{i,2+j}) max(lanes{i,2+j}) ];
        ds = scale(2)-scale(1);
        scale(1) = scale(1)-0.1*ds; scale(2) = scale(2)+0.1*ds; 

        
        set(gca, 'XLim', scale)
        legend(channel{j})
    end

    
    print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_lane_' num2str(i) '.png']); %save figure

    
end


























    
    
    
    
    
    
