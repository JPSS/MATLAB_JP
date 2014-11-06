%%
clc
clear all 
close all

path0 = cd; addpath(path0); display(['Using path: ' path0 ])
matlab_dir = userpath; matlab_dir = matlab_dir(1:end-1);
run([matlab_dir filesep 'my_prefs.m'])
scrsz = get(0,'ScreenSize');
    cur_fig = figure('Visible','on','OuterPosition',[ 1 scrsz(4) scrsz(4)*0.5 scrsz(4)*0.5], 'PaperPositionMode', 'auto'); % figure('Visible','off');%left bottom width height


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
    angle = tmp(:, 6);
    
    s_tmp{1,1} = angle;

    name = inputdlg({'Length of spacer:'}, 'Length of spacer' , 1, {'10'} );
    s_tmp{1, 2} = str2double(name{1});
    
    

    


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
    
    print(cur_fig, '-djpeg' , '-r300', [pname fname '_plot.jpg']); %save figure

    structure = [s_tmp; structure];
    button = questdlg('More data?','More data','More','Enough', 'Enough');
    more = strcmp(button,'More');    
    
    
end


%% PLOT DISTRIBUTIONS
xhist_angle = 0:1:35;
xhist_rad = sin(xhist_angle(1)*pi/180/2):0.01:sin(xhist_angle(end)*pi/180/2);
for i=1:size(structure, 1)
   cur_fig = figure('Visible','on','OuterPosition',[ 1 scrsz(4) scrsz(4)*0.6 scrsz(4)*0.8], 'PaperPositionMode', 'auto'); % figure('Visible','off');%left bottom width height

    %angle
    subplot(2, 1, 1)
    [p n] = hist_fit(structure{i,1}, xhist_angle);
    bar(xhist_angle, n, 'b'), hold on
    plot(xhist_angle, gauss1d(p, xhist_angle), 'r')
    xlabel('Angle', 'Fontsize', 12)
    ylabel('Frequency','Fontsize', 12)
    title({['Angle-Frequency Distribution of 4sp1-' num2str(structure{i,2})], ['mean = ' num2str(mean(structure{i,1})) ' +- ' num2str(std(structure{i,1})/sqrt(length(structure{i,1})))]}, 'Fontsize', 14)
    structure{i,3} = p;
    
    %angle
    subplot(2, 1, 2)
    rad = sin(structure{i,1}*pi/180/2);
    [p n] = hist_fit(rad, xhist_rad);
    bar(xhist_rad, n, 'b'), hold on
    plot(xhist_rad, gauss1d(p, xhist_rad), 'r')
    xlabel('sin(\alpha/2)', 'Fontsize', 12)
    ylabel('Frequency','Fontsize', 12)
    title({['Angle-Frequency Distribution of 4sp1-' num2str(structure{i,2})], ['mean = ' num2str(mean(rad)) ' +- ' num2str(std(rad)/sqrt(length(rad)))]}, 'Fontsize', 14)
    structure{i,4} = p;

end

%% PLOT distributions AND length vs. sin(alpha/2)
my_colors = [1 0 0; 0 1 0; 0 0 1];
for i=1:3
    subplot(1, 5, i)
    [p n] = hist_fit(sin(structure{i,1}*pi/180/2), xhist_rad);
    barh(xhist_rad, n, 'FaceColor', my_colors(i,:)), hold on
    plot( gauss1d(p, xhist_rad), xhist_rad,'k')
    set(gca, 'YLim', [xhist_rad(1) xhist_rad(end)])
    if i ==1
        ylabel('sin(\alpha/2)', 'Fontsize', 12)
    else
        set(gca, 'YTick', [])
    end
    xlabel('Frequency','Fontsize', 12)
    title({[num2str(structure{i,2}) 'bp']}, 'Fontsize', 14)
    legend({[ 'data(' num2str(sum(n)) ')' ], 'fit'}, 'Fontsize', 12)

    
    
    
    mean_rad(i) = mean(sin(structure{i,1}*pi/180/2));
    error_rad(i) = std(sin(structure{i,1}*pi/180/2)) / sqrt(length(sin(structure{i,1}*pi/180/2)));
    n_bp(i) = structure{i,2};
    
    subplot(1, 5, 4:5)
    errorbar(n_bp(i), mean_rad(i), error_rad(i), '.','Color', my_colors(i,:), 'MarkerSize', 25), hold on
    
end


p = polyfit(n_bp, mean_rad, 1);

subplot(1, 5, 4:5)
%errorbar(n_bp, mean_rad, error_rad, '.b'), hold on
plot([0:20], p(1)*[0:20]+p(2), 'k--'), hold on
set(gca, 'YLim', [xhist_rad(1) xhist_rad(end)])
set(gca, 'YAxislocation', 'right', 'Fontsize', 12)
legend({'10bp mean +- error', '14bp mean +- error', '18bp mean +- error','fit'}, 'Fontsize', 12)
xlabel('Length of spacer [bp]', 'Fontsize', 12)
ylabel('sin(\alpha/2)', 'Fontsize', 12)

%%
x = zeros(size(structure, 1), 1) ;
alpha =zeros(size(structure, 1), 2) ;

alpha_fit =zeros(size(structure, 1), 2) ;
for i=1:size(structure,1)
    x(i) = structure{i,2};
    alpha(i,1) = mean(structure{i,1});
    alpha(i,2) = std(structure{i,1}) / sqrt(length(structure{i,1}));
    
    alpha_fit(i,1) = structure{i,3}(1);
    alpha_fit(i,2) = structure{i,3}(2);
end
s_alpha(:,1) = sin(alpha(:,1)*pi/180/2);
s_alpha_fit(:,1) = sin(alpha_fit(:,1)*pi/180/2);

s_alpha(:,2) = (alpha(:,2)/2).*(pi./180).*cos(alpha(:,1)*pi/180/2);
s_alpha_fit(:,2) = (alpha_fit(:,2)/2).*(pi./180).*cos(alpha_fit(:,1)*pi/180/2);


p = polyfit(x, s_alpha(:,1), 1);

p_fit = polyfit(x, s_alpha_fit(:,1), 1);





%%
subplot(1, 5, 4:5)
errorbar(x, s_alpha(:,1), s_alpha(:,2), '.b'), hold on
plot(x, p(1)*x+p(2), 'b--'), hold on
errorbar(x, s_alpha_fit(:,1), s_alpha_fit(:,2), '.r'), hold on
plot(x, p_fit(1)*x+p_fit(2), 'r--')

legend('Data', 'fit' , 'data from fit', 'fit')
xlabel('Length of spacer [bp]', 'Fontsize', 12)
ylabel('sin(\alpha/2)', 'Fontsize', 12)

%set(gca , 'YLim', [0 0.4]);
set(gca, 'YLim', [0 0.3])
nbins = 10;

xhist = 0: 0.01 :0.4;


subplot(1, 5, 1)
%cs = varycolor(size(structure, 1));

a = sin(structure{3,1}*pi/180/2);
n = hist(a, xhist);
barh(xhist, n,  'facecolor','g' ), hold on
plot( sin(gauss1d(structure{3,3}, xhist)*pi/180/2), xhist, 'k')
set(gca, 'YLim', [0 0.3])

subplot(1, 5, 2)

a = sin(structure{2,1}*pi/180/2);
n = hist(a, xhist);
barh(xhist, n,  'facecolor','r' ), hold on
p = [structure{2,3}(1:2)*pi/180/2 structure{2,3}(3) ]
plot( gauss1d(p, xhist), xhist, 'k')
set(gca, 'YLim', [0 0.3])

subplot(1, 5, 3)

a = sin(structure{1,1}*pi/180/2);
n = hist(a, xhist);
barh(xhist, n,  'facecolor','b' ), hold on
p = [structure{1,3}(1:2)*pi/180/2 structure{1,3}(3) ]
plot( gauss1d(p, xhist), xhist, 'k')
set(gca, 'YLim', [0 0.3])


%%
for i=1:size(structure, 1)
    a = sin(structure{i,1}*pi/180/2);
    
    n = hist(a, xhist);
    barh(xhist, n,  'facecolor',cs(:,i) , 'edgecolor',cs(:,i)), hold on
   
end
set(gca , 'YLim', [0 0.4]);


%%
print(cur_fig, '-djpeg' , '-r300', [pname fname '_all.jpg']); %save figure





save([pname fname '_fit.txt'], 'p' ,'-ascii')

cd(path0)


