%%
clear all, close all, clc
run('my_prefs')
path0 = cd;

%% load image
cd(data_dir)
[filename pathname]=uigetfile('*.tif','Select image:');
cd(pathname)
cd(path0)


%% create output folder
pname = inputdlg({'Output folder and prefix:'}, 'Output folder and prefix' , 1, {filename(1:size(filename,2)-4)} );
prefix_out = pname{1};
path_out = [pathname prefix_out ];
mkdir(path_out);
%% load and bg correct images
img = double(imread([pathname filesep filename]));  %load
plot_image_ui(img)
button = questdlg('Rotate?','Rotate','Rotate','No','No');
if strcmp(button,'Rotate') %load old data
    img = imrotate(img, -90);
end
close all
img_bg = bg_correct_ui(img, 'Background correction');    %bg correct
close all

%% select bands by hand
imagesc(img_bg), colormap gray
options.WindowStyle='normal';
prompt={'How many bands'};
def={'3'};
tmp = inputdlg(prompt, 'How many bands', 1, def, options);
n_bands = str2double(tmp(1));
close all
%%
[I, areas] = integrate_areas({img_bg}, n_bands, 1)

%%
%colormap('gray'), colorbar

%imagesc(img_bg)
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual');
colormap('gray'), colorbar

imagesc(img_bg)
%plot_image_ui(img_bg)
for i=1:n_bands
   rectangle('Position', areas(i,:), 'EdgeColor', 'r')
   text(double(areas(i,1)), double(areas(i,2)), num2str(i), 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'top')
end 
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_bands.eps']); %save figure


%%
save([path_out filesep prefix_out '_data'] )

%% PLOTTING STUFF

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I)
%set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',names,'Fontsize', 12)
%set(gca, 'XLim', [0 size(lanes,1)+1])
%xticklabel_rotate([1:size(lanes,1)],90,names)

set(gca, 'XTick', 1:n_bands)
xlabel('Band')
ylabel('Intensity')
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_IntesityRaw.eps']); %save figure


%%

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(  (I(3:3:end)+I(2:3:end) ) ./ ( I(1:3:end)+I(3:3:end)+I(2:3:end) ))
%set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',names,'Fontsize', 12)
%set(gca, 'XLim', [0 size(lanes,1)+1])
%xticklabel_rotate([1:size(lanes,1)],90,names)

set(gca, 'XTick', 1:n_bands)
xlabel('Lane')
ylabel('norm. Intensity 1-fach + 2-fach / alles')
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_labeled_vs_all.eps']); %save figure



%%

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(  (I(2:3:end) ) ./ ( I(1:3:end)+I(2:3:end) ))
%set(gca, 'XTick',1:size(lanes,1), 'XTickLabel',names,'Fontsize', 12)
%set(gca, 'XLim', [0 size(lanes,1)+1])
%xticklabel_rotate([1:size(lanes,1)],90,names)

set(gca, 'XTick', 1:n_bands)
xlabel('Lane')
ylabel('norm. Intensity 1-fach  / 0-fach+1-fach')
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_1-fach-labeled_vs_0+1.eps']); %save figure






    
%%
close all
c = [2  1.5 1 0.75 0.5 0.25 0.1]'*1000; % ug/ml
V = ones(size(c)); 
n = c .* V;
[c] = polyfit(n, I(1:7), 1);
plot(n, I(1:7), 'b.', 'MarkerSize', 15), hold on
plot(0:n(1), c(2)+c(1)*[0:n(1)], 'r--')
hline(I(8:end));

%%
m_BSA = 66.5;
m_tetR = 51380; %Da, g/mol
m_K271 = 78.29e3; %Da, g/mol
m_v3 = 51.23e3;

conc = ((I(8:end)-c(2))/c(1)).*[1 10 20 1 10 20]' / 1000 %g/L


(conc(1:3)/(m_K271))*1e6%uMol/l
(conc(4:end)/(m_v3))*1e6%uMol/l

%%
semilogy(n, I, 'b.', 'MarkerSize', 15)
%%
dlmwrite([path_out filesep 'data.txt'], [c V I], '\t')


%%
save([path_out filesep prefix_out '_data'] )
display('done')

%%
figure(2)
bar(I)
    




%%


close all
c = [1000 500 250]';
V = (10/12)*[11 11 11]';
n = c .* V;
[cf, c_err, cov] = fit_linear(n, I(1:3));
plot(n, I(1:3), 'b.', 'MarkerSize', 15), hold on
plot(n(end):n(1), cf(2)+cf(1)*[n(end):n(1)], 'r--')

%%
figure(2)
plot( [500 1000 700], I)

    
    