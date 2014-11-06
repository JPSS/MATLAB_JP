function [ param, param_err ] = detect_peak(x, y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% select area of interest
close all
plot(x,y)
xlim = [x(1) x(end)];
set(gca, 'XLim', xlim)
ylim = get(gca, 'Ylim');

h = imrect;
pos = int32(wait(h)); % [xmin ymin width height]

x1 = pos(1);
x2 = pos(1)+pos(3);
close all

%% fit a gaussian to selected area 
sigma = 5;

dx = x(2)-x(1);
i1 = double(round(   (x1-x(1))/dx  )+1);
i2 = double(round(   (x2-x(1))/dx  )+1);
[y_max x_max] = max(y(i1:i2));
x_max = x_max + pos(1);
p_init = double([x_max sigma y_max ]);

options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',
[param,r,J,COVB,mse] = nlinfit(x(i1:i2),  y(i1:i2) ,@gauss1d,p_init, options);
ci = nlparci(param,r,'Jacobian',J);
param_err = abs( ci(:,1)' - param);


x_plot = x(1):(x(2)-x(1))/10:x(end);
plot(x, y, 'g', x_plot,  gauss1d(param, x_plot), 'r--' )


param(2) = abs(param(2));
end

