%% startup
clc
clear all 
close all

run('my_prefs.m')
path0 = cd; 

%% get file
[fname pname] = uigetfile('.txt', 'Select txt-file with angles');
cd(path0)

%% get data from file
[data, raw_data] = get_structure([pname fname], 'GFP 3-132', 10);

%%
dlmwrite([pname 'mean_std_err.txt'], [mean(data{1}) std(data{1}) std(data{1})/sqrt(length(data{1}))], 'delimiter', '\t')
dlmwrite([pname 'slice_angle.txt'], [data{2} data{1}], 'delimiter', '\t')


%% Plot distribution
x_min = floor(min(data{1})/10)*10;
x_max = ceil(max(data{1})/10)*10;
dx = 10;
xhist = x_min:dx:x_max;
xplot = x_min:dx/100:x_max;
n  = hist(data{1}, xhist);

close all
fig_dim =[10 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

h = bar(xhist, n);
set(gca, 'YLim', [0 max(n)*1.5])
legend(h, {'Data'})
set(gca, 'XLim', [x_min x_max])
xlabel('Angle [deg]'), ylabel('Frequency')

print(cur_fig, '-dtiff', '-r600' , [pname 'Histogram_01_bin=' num2str(dx) '.tif'])

%% Plot with ONE gauss
close all
fig_dim =[10 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(2,1, 1)
h(1) = bar(xhist, n./sum(n)), hold on;
h(2)= plot(xplot, gauss1d([mean(data{1}) std(data{1}) dx/(sqrt(2*pi)*std(data{1}))], xplot), 'r')
set(gca, 'YLim', [0 max(n./sum(n))*1.5])
legend(h, {'Data', 'Gaussian-Distribution'})
set(gca, 'XLim', [x_min x_max])
 ylabel('rel. Frequency')

subplot(2,1,2)
plot(xhist, gauss1d([mean(data{1}) std(data{1}) 1/(sqrt(2*pi)*std(data{1}))], xhist) - n./sum(n), 'b.-')
%set(gca, 'YLim', [0 max(n./sum(n))*1.5])
set(gca, 'XLim', [x_min x_max])
xlabel('Angle [deg]'), ylabel('Residuals')

print(cur_fig, '-dtiff', '-r600' , [pname 'Histogram_one-gauss_bin=' num2str(dx) '.tif'])


%% plot energy function
x_min = 2*47.5*0.34*sin(floor(min(data{1})/10)*10*pi/180/2);
x_max = 2*47.5*0.34*sin(ceil(max(data{1})/10)*10*pi/180/2);
dx = 2*47.5*0.34*sin(3*pi/180/2);
%xhist = [x_min:dx:x_max];
xplot = [x_min:dx/100:x_max];
m  = hist(2*47.5*0.34*sin(data{1}*pi/180/2), xhist);
n = m ./ sum(m);

close all
fig_dim = [20, 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
subplot(1,3,1)
fig_dim =[10 10];
h = bar(xhist, n);
%set(gca, 'YLim', [0 max(n)*1.5])
legend(h, {'Data'})
set(gca, 'XLim', [x_min x_max])
xlabel('Length [bp]'), ylabel('rel. Frequency')

subplot(1,3,2)

%E = -log(n);
i0 =2;
i1 = length(E)-2;
plot(xhist, E, 'b.'), hold on
plot(xhist(i0:i1),E(i0:i1), 'r.')
p = polyfit(xhist(i0:i1),E(i0:i1), 3 );
plot(xplot, polyval(p, xplot) , 'g-') %p(3)+p(2)*xplot.^1+p(1)*xplot.^2
set(gca, 'YLIM', [1 7])

subplot(1,3,3)
dxplot = (xplot(2)-xplot(1))/2;
plot(xplot(2:end)+dxplot/2, diff(kT*polyval(p, xplot))/dxplot , 'g-') %p(3)+p(2)*xplot.^1+p(1)*xplot.^2


%%
x0 = p(2)/p(1)/(-2)
offset = p(3)-x0^2*p(1)
a = p(1)

subplot(1,3,3)
%%
close all
plot(xhist(2:end)+dx/2 , diff(kT*E)/dx, 'r',xhist, E, 'b' ) 

%%
subplot(1,3,3)
kB = 1.3806488e-23; % J/K
T = 297; %K
kT = kB*T*1e21; % pN*nm
plot(xplot/0.34, 2*kT*a*(xplot-x0))
xlabel('Lenght of spacer [bp]')
ylabel('Force [pN]')
%% Two gauss
%{
xhist = 0:2:120;
xplot = 0:0.1:120;
n  = hist(data{1}, xhist);


options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',

p_init = [49 4 5 60 4 55];

[p,r,J,COVB,mse] = nlinfit(xhist, n,@two_gauss1d, p_init, options);
ci = nlparci(p,r,'Jacobian',J);

p(2) = abs(p(2));
p(5) = abs(p(5));
p_err =  abs(p-ci(:,1)')  ;

p
%
close all
fig_dim =[10 7.5];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(xhist, n), hold on
plot(xplot, two_gauss1d(p , xplot), 'r', 'Linewidth', 3)
set(gca, 'XLim', [30 100])
xlabel('Angle [deg]'), ylabel('Frequency')

print(cur_fig, '-dtiff', '-r500' , [pname 'Histogram_02.tif'])

close all
fig_dim =[10 7.5];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(xhist, n), hold on
plot(xplot, gauss1d(p(1:3) , xplot), 'r', xplot, gauss1d(p(4:6) , xplot), 'g','Linewidth', 3)
set(gca, 'XLim', [30 100])
xlabel('Angle [deg]'), ylabel('Frequency')

print(cur_fig, '-dtiff', '-r500' , [pname 'Histogram_03.tif'])
%}
%% One gauss old
%{ 
FIT ONE GAUSS
dx = 2;
xhist = 0:dx:120;
xplot = 0:0.1:120;
n  = hist(data{1}, xhist);


options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',

p_init = [mean(data{1}) std(data{1}) max(n)];

[p,r,J,COVB,mse] = nlinfit(xhist, n,@gauss1d, p_init, options);
ci = nlparci(p,r,'Jacobian',J);

p(2) = abs(p(2));
p_err =  abs(p-ci(:,1)')  ;

p
p_err

    
%
close all
fig_dim =[10 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);



subplot(4, 1, 1:3)
h = zeros(4,1);
h(1) = bar(xhist, n); hold on;
set(gca, 'YLim', [0 max(n)*1.5])
h(2) = plot(xplot, gauss1d(p , xplot), 'r','Linewidth', 2); hold on;

%plot(xplot, gauss1d(p , xplot), 'r',xplot, dx*gauss1d([mean(data{1}) std(data{1}) 1/sqrt(2*pi*std(data{1}).^2)] , xplot), 'g','Linewidth', 2); hold on;
h(3) = vline(p(1), 'r-');
h(4) = vline(mean(data{1}), 'r--');
legend(h,{'Data', 'Gaussian fit', 'Mean of fit', 'Mean of data'})
set(gca, 'XLim', [0 60])
ylabel('Frequency')


subplot(4, 1, 4)
plot(xhist, n-gauss1d(p , xhist), 'b.-'), hold on
hline(0, 'k')
set(gca, 'XLim', [0 40])
xlabel('Angle [deg]'), ylabel('Frequency')
ylabel('Residuals')

print(cur_fig, '-dtiff', '-r600' , [pname 'Histogram_01_bin=' num2str(dx) '.tif'])
%}

%% write angles to file


%% wirte images to output file
cd(pname)
[fname2 pname2] = uigetfile('.png', 'Select img-file');
cd(path0)
prefix = fname2(1:end-7)
%%
n_img = size(raw_data,1)/2;
close all
colormap gray
arrows = zeros(n_img*2,4);
for i=1:n_img*2
   
   if 0 <= raw_data(i,6) && raw_data(i,6) < 90
       x = raw_data(i,2);
       y = raw_data(i,3)+raw_data(i,5); %oben links + height
   else
       if 90 <= raw_data(i,6) && raw_data(i,6) < 180
           x = raw_data(i,2)+raw_data(i,4);
           y = raw_data(i,3)+raw_data(i,5); 
       else
           if -180 <= raw_data(i,6) && raw_data(i,6) < -90
               x = raw_data(i,2)+raw_data(i,4);
               y = raw_data(i,3); 
           else
               if -90 <= raw_data(i,6) && raw_data(i,6) < 0
                   x = raw_data(i,2);
                   y = raw_data(i,3); 
               end
           end
       end
   end   
   r = sqrt(raw_data(i,4).^2+raw_data(i,5).^2);
   dx = r*cos(raw_data(i,6)*pi/180);
   dy = -r*sin(raw_data(i,6)*pi/180);
   arrows(i,:) = [x+1 y+1 dx dy];
end

%%
fig_dim =[10 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

for i=0:10%n_img-1
   %load image
   img = imread([pname2 prefix sprintf('.%02i.tif', raw_data(2*i+1,7)-1) ]);
   img2 = max(img(:))-img;

   imagesc(img2), colorbar, hold on
   axis image
   colormap gray

   quiver(arrows(2*i+1,1), arrows(2*i+1,2), arrows(2*i+1,3), arrows(2*i+1,4), 0, 'r', 'linewidth', 1.5)
    quiver(arrows(2*i+2,1), arrows(2*i+2,2), arrows(2*i+2,3), arrows(2*i+2,4), 0, 'r', 'linewidth', 1.5)

    
    hold off
    
    
    print(cur_fig, '-dpng','-r300' , [pname2 prefix sprintf('_arrow.%02i.tif', raw_data(2*i+1,7)-1) ]); %save figure

    % pause

end
    

%%
close all
for i=0:0%n_img-1
   %load image
   img = imread([pname2 prefix sprintf('.%02i.tif', raw_data(2*i+1,7)-1) ]);
   img2 = max(img(:))-img;
   %hist(double(img(:)), 100), hold on
   imagesc(img2, [min(double(img2(:))) 2*max(double(img2(:)))]), colorbar, hold on
      colormap gray


   
   
end
%%









