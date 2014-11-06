
%%
clear all
close all
clc
run('my_prefs')


%% for sides form 

x_a=7; %1st peak position
y_a=5; %1st peak position
s=1.3;  %sigma 
a=10; %amplitude 
d=4.8; %distance

x_b=x_a;
y_b=y_a+d;

N=14; % box dimension (pxl)
z = zeros(N,N);


 for i=1:N % i is 1pxl step
     for j=1:N % j is 1pxl step
     
 z(i,j) = a.*exp(-((((i-x_a)^2)+(j-y_a)^2))/(2*(s)^2)) + 2*a.*exp(-((((i-x_b)^2)+(j-y_b)^2))/(2*(s)^2));
  
     end
 end

 
%figure
surf(z); %plot matrix in 3d, one color
%mesh(z); %plot matrix in 3d with colour gradient in z


scrsz = get(0,'ScreenSize');
fig_dim =[N N];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

imagesc(z)

colormap gray
%colormap gray / jet
colorbar
axis image
title( ['d='  num2str(d) 'px sigma=' num2str(s) 'px' ])

path_out = '/Users/letiziameregalli/Documents/PSF_simulation/';
print(cur_fig, '-dtiff', '-r300',   [path_out  'sides_d='  num2str(d) 'px_sigma=' num2str(s) 'px_gray.tif' ])
