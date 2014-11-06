%%
clc
clear all 
close all

path0 = cd; addpath(path0); display(['Using path: ' path0 ])
matlab_dir = userpath; matlab_dir = matlab_dir(1:end-1);
run([matlab_dir filesep 'my_prefs.m'])
scrsz = get(0,'ScreenSize');
   % cur_fig = figure('Visible','on','OuterPosition',[ 1 scrsz(4) scrsz(4)*0.5 scrsz(4)*0.5], 'PaperPositionMode', 'auto'); % figure('Visible','off');%left bottom width height
%%
structure = cell(0,4);
more = 1;
while more
    
   
    s_tmp = cell(1,4);
    [fname pname] = uigetfile('.txt', 'Select txt-file with angles');
    cd(pname)
    cd ..
    cd ..
    cd ..
    tmp = load([pname fname]);
    
    
    %%
    img_path = pname;
    img_names = cell(0,1);
    
    dangle = zeros(size(tmp,1)/2, 1);
    for i=1:2:size(tmp,1)
        img_names = [img_names ; ['sideview.' sprintf('%02i', tmp(i,7)-1) '.tif']];
        
        
        dangle((i+1)/2) = abs(tmp(i,6) - tmp(i+1,6)) ;
        if dangle((i+1)/2) > 180
            dangle((i+1)/2) = 360 - dangle((i+1)/2);
        end
    end
    %%
    
    s_tmp{1,1} = dangle;
    s_tmp{1, 2} = img_names;
    s_tmp{1, 3} = img_path;
    
    
    
    name = inputdlg({'Name of stack:'}, 'Name' , 1, {'dsLoop'} );
    s_tmp{1, 4} = name{1};
    
    

    

    %{
    nbins = max(20, sqrt(length(angle)));
    xhist = 0:(max(angle)-0)/(nbins-1):max(angle);

    [p n] = hist_fit(angle, xhist);

    s_tmp{1, 3} = p;
  
    
    
    bar(xhist, n, 'b'), hold on
    plot(xhist, gauss1d(p, xhist), 'r')
    legend([num2str(length(angle)) ' data points'], ['fit max at ' num2str(p(1))])
    xlabel('Angle', 'Fontsize', 12)
    ylabel('Frequency','Fontsize', 12)
    title({['Angle-Frequency Distribution of 4sp1-' num2str(s_tmp{1,2})], ['mean = ' num2str(mean(angle)) ' +- ' num2str(std(angle)/sqrt(length(angle)))]}, 'Fontsize', 14)
    hold off
    
    %print(cur_fig, '-djpeg' , '-r300', [pname fname '_plot.jpg']); %save figure
    %}
    
    structure = [s_tmp; structure];
    button = questdlg('More data?','More data','More','Enough', 'Enough');
    more = strcmp(button,'More');    
    
    
end

%%
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 13.7 20]); 

%xlabel('Wavelength [nm]', 'FontSize', 12,'FontName','Times New Roman')
%ylabel('Molar exctinction coefficient [M^{-1}cm^{-1}]', 'FontSize', 12,'FontName','Times New Roman' )





xhist = 0:2:15;
n = zeros(size(structure,1), length(xhist));
for i=1:size(structure, 1)
    n(i, :) = hist(structure{i,1}, xhist);
    
    
    %plot histogram
    subplot(size(structure, 1), 1, i)
    bar(xhist, n(i, :))
    title(structure{i,4},'FontSize', 14)
    xlabel('Angle [deg]','FontSize', 12)
    ylabel('Frequency','FontSize', 12)
    legend([num2str(length(structure{i,1})) ' data points'],'FontSize', 12)
           set(gca, 'box', 'off')

    %plot energy landscape
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none',  'XColor','r','YColor','r');
    set(ax2, 'XTick', get(ax1, 'XTick'))
    set(ax2, 'XLim', get(ax1, 'XLim'))
    %set(ax2, 'YTick', get(ax1, 'YTick'))

       set(ax2, 'XTickLabel', [])
       set(ax2, 'XTick', [])
       set(ax2, 'XColor', 'w')
    line(xhist, -log(n(i,:)/sum(n(i,:)) ), 'Color','r', 'LineWidth', 2)% ,'Parent',ax2)
    ylabel('E [k_BT]','FontSize', 12)
    
end
%print(cur_fig, '-depsc','-cmyk','-loose' , 'test'); %save figure
print(cur_fig, '-depsc','-loose' , 'test'); %save figure
%%

figure(2)
xhist = 2:0.1:8;
n = zeros(size(structure,1), length(xhist));
for i=1:size(structure, 1)
    n(i, :) = hist(sqrt(2)*3.73*sqrt(1+cos(structure{i,1}*pi/180)), xhist);   % length of spring = ... [nm]
    
    
    %plot histogram
    subplot(size(structure, 1), 1, i)
    bar(xhist, n(i, :))
    title(structure{i,4})
    xlabel('Extension of spring [nm]')
    ylabel('Frequency')
    legend([num2str(length(structure{i,1})) ' data points'])
    
    %plot energy landscape
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none',  'XColor','r','YColor','r');
    set(ax2, 'XTick', get(ax1, 'XTick'))
    set(ax2, 'XLim', get(ax1, 'XLim'))
    %set(ax2, 'YTick', get(ax1, 'YTick'))

    line(xhist, -log(n(i,:)/sum(n(i,:)) ), 'Color','r')% ,'Parent',ax2)
    ylabel('E [k_BT]')
    
end




figure(3)
xhist = 1:6:80;
n = zeros(size(structure,1), length(xhist));
for i=1:size(structure, 1)
    n(i, :) = hist(2*47*sin(structure{i,1}*pi/180/2), xhist);   % length of spring = ... [nm]
    
    
    %plot histogram
    subplot(size(structure, 1), 1, i)
    bar(xhist, n(i, :))
    title(structure{i,4})
    xlabel('Extension of spacer/sample [bp]')
    ylabel('Frequency')
    legend([num2str(length(structure{i,1})) ' data points'])
    
    %plot energy landscape
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none',  'XColor','r','YColor','r');
    set(ax2, 'XTick', get(ax1, 'XTick'))
    set(ax2, 'XLim', get(ax1, 'XLim'))
    %set(ax2, 'YTick', get(ax1, 'YTick'))

    line(xhist, -log(n(i,:)/sum(n(i,:)) ), 'Color','r')% ,'Parent',ax2)
    ylabel('E [k_BT]')
    
    
    
    
end



%{
subplot(2, 1, 2)
bar(xhist, n(2, :))
title(structure{2,4})
xlabel('Angle [deg]')
ylabel('Frequency')
legend([num2str(length(structure{1,1})) 'data points'])
%}



%%
