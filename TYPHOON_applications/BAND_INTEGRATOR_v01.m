%%
clear all, close all, clc
run('my_prefs')
path0 = cd;

%% load image
nimg_input = inputdlg({'Number of images:'}, 'Number of images' , 1, {'1'} );
n_img = str2double(nimg_input{1});

filenames = cell(n_img, 1);
pathnames = cell(n_img, 1);

last_dir = data_dir;
for i=1:n_img
    cd(last_dir)
    [filenames{i} pathnames{i}]=uigetfile('*.tif','Select image:');
    last_dir = pathnames{i};
end
cd(path0)

%% create output folder
pname = inputdlg({'Output folder and prefix:'}, 'Output folder and prefix' , 1, {filenames{1}(1:size(filenames{1},2)-4)} );
prefix_out = pname{1};
path_out = [pathnames{1} prefix_out ];
mkdir(path_out);

%% load and bg correct images
images = cell(n_img, 1);
img_bg = cell(n_img, 1);
for i=1:n_img
    images{i} = double(imread([pathnames{i} filesep filenames{i}]));  %load
    plot_image_ui(images{i});
    button = questdlg('Rotate?','Rotate','Rotate','No','No');
    if strcmp(button,'Rotate') %load old data
        images{i} = imrotate(images{i}, -90);
    end
    close all
    img_bg{i} = bg_correct_ui(images{i}, 'Background correction');    %bg correct
    close all    
end

%% select bands by hand
imagesc(img_bg{1}), colormap gray
options.WindowStyle='normal';
prompt={'How many bands'};
def={'3'};
tmp = inputdlg(prompt, 'How many bands', 1, def, options);
n_bands = str2double(tmp(1));
close all

%% integrate areas, all areas have the same size
[I, areas] = integrate_areas(img_bg, n_bands, 1); %cell of images, number of bands, 1=all bands habe the same size

%% assining names
names = cell(n_bands,1);
for i=1:n_bands
    band_name = inputdlg({'Name of band:'}, 'Band names' , 1, {['Band ' num2str(i)]} );
    names{i} = band_name{1};
    
end


%%
%colormap('gray'), colorbar

%imagesc(img_bg)
for i=1:n_img
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual');
    colormap('gray'), colorbar

    imagesc(img_bg{i}), colorbar, axis image, colormap gray,
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

%%
save([path_out filesep prefix_out '_data'] )

dlmwrite([path_out filesep 'data.txt' ], I, 'delimiter', '\t')


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

bar(1:n_bands, I)
set(gca, 'XTick',1:n_bands, 'XTickLabel',names,'Fontsize', 12)
xticklabel_rotate([1:n_bands],90,names)
ylabel('Intensity')
%legend({'GFP', 'DNA'})
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_IntesityRaw.png']); %save figure
print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_IntesityRaw.tif']); %save figure

%% PLOTTING STUFF
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I(1:2:n_bands,:)+I(2:2:n_bands,:))

set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names(2:2:n_bands,1),'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names(2:2:n_bands,1))

set(gca, 'XTick', 1:n_bands)
%ylabel('Sum of labeled- and unlabeled-band')
ylabel('Sum of monomer and dimer band')
legend({'GFP', 'DNA'})

print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_sum.png']); %save figure
print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_sum.tif']); %save figure

%%
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


bar(I(2:2:n_bands,:)./(I(1:2:n_bands,:)+I(2:2:n_bands,:)))
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names(2:2:n_bands,1),'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names(2:2:n_bands,1))

set(gca, 'XTick', 1:n_bands/2)
set(gca, 'YLim', [0 1])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('[dimer] / ([dimer] + [monomer])')
%ylabel('[GFP] / ([GFP] + [GFP+DNA])')

%legend({'GFP-channel', 'DNA-channel'})

print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_normalizedIntensities.png']); %save figure
print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_normalizedIntensities.tif']); %save figure


%%
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


bar(I(2:2:n_bands,1)./(I(1:2:n_bands,1)+I(2:2:n_bands,1)))
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names(2:2:n_bands,1),'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names(2:2:n_bands,1))

set(gca, 'XTick', 1:n_bands/2)
set(gca, 'YLim', [0 1.5])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('[GFP+DNA] / ([GFP] + [GFP+DNA])')

%legend({'GFP-channel', 'DNA-channel'})
legend({'GFP-channel'})


print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_normalizedIntensities2.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_normalizedIntensities2.png']); %save figure


%% ratio of bands
close all
fig_dim =1.5*[10 7.5];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I(:,2)./I(:,1))
set(gca, 'XTick',1:n_bands, 'XTickLabel',names,'Fontsize', 10)
xticklabel_rotate([1:n_bands],90,names)

set(gca, 'XTick', 1:n_bands)


for i=1:n_bands
    text(i,0.5, [num2str(round(I(i,2)./I(i,1)*10)/10)], 'Color', 'white', 'FontSize', 15, 'HorizontalAlignment','center')
end


set(gca, 'YLim', [0 1.5])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('Ratio of [Reporter] / [GFP]')

%legend({'GFP-channel', 'DNA-channel'})
%legend({'GFP-channel'})


print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_relative.eps']); %save figure
print(cur_fig, '-dpng2' , [path_out filesep prefix_out '_relative.png']); %save figure


%% ratio of bands
close all
fig_dim =2.5*[10 7.5];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I(2:2:end,2)./I(1:2:end,1))
Iset(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 10)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands/2)


for i=1:2:n_bands
    text(1+(i-1)/2,0.5, [num2str(round(I(i+1,2)./I(i,1)*100)/100)], 'Color', 'red', 'FontSize', 15, 'HorizontalAlignment','center')
end


%set(gca, 'YLim', [0 1.5])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('Ratio of [DNA] / [GFP]')

%legend({'GFP-channel', 'DNA-channel'})
%legend({'GFP-channel'})


print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_relative.eps']); %save figure
print(cur_fig, '-dpng', '-loose'  , [path_out filesep prefix_out '_relative.png']); %save figure
print(cur_fig, '-dtiff','-r300' , '-loose' , [path_out filesep prefix_out '_relative.tif']); %save figure







%% Old stuff
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

c = [5 7.5 10 15 20]'; % uM
I_standard = I(7:11);

[cf, c_err, cov] = fit_linear(c, I_standard);
plot(c, I_standard, '.', 'MarkerSize', 15), hold on
plot(0:c(end), cf(2)+cf(1)*[0:c(end)], 'r--')

(I-cf(2))/cf(1)




%%
for i=1:12
    if i<4
        I_norm(i) = sum(I(1:3));
    else
        if i <7
            I_norm(i) = sum(I(4:6));
        else
            if i <10
                I_norm(i) = sum(I(7:9))
            else
                I_norm(i) = sum(I(10:12))
            end
        end
    end
end

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);



bar(I(1:12) ./ I_norm' )
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_relativeIntesity.eps']); %save figure

    


%% 
close all
fig_dim =1.5*[10 7.5];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

plot(1:16, I(1,:)./I(1,1), 'r.-', 1:16, I(2,:)./I(2,1), 'b.-', 'Markersize', 15)
legend(names(1:2), 'location', 'southeast')
set(gca, 'YLim', [0 1.5])
xlabel('Scan nr.')
ylabel('Normalized intensity')

print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_trace.eps']); %save figure
print(cur_fig, '-dpng', '-loose'  , [path_out filesep prefix_out '_trace.png']); %save figure
print(cur_fig, '-dtiff','-r500' , '-loose' , [path_out filesep prefix_out '_trace.tif']); %save figure



%%


I = zeros(n_bands, n_img);
for i=1:n_bands
    pos = areas(i,:);
    for j = 1:n_img
        I(i,j) = sum(sum((img_bg{j}(pos(2):pos(2)+pos(4)  , pos(1):pos(1)+pos(3))))); %integrate
    end
    close all
end
    