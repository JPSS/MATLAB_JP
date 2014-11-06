function [ p_fast, p_slow, p_fast_err, p_slow_err] = fit_2peaks(x, y )
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
sigma = 8;

dx = x(2)-x(1);
i1 = double(round(   (x1-x(1))/dx  )+1);
i2 = double(round(   (x2-x(1))/dx  )+1);
[y_max x_max] = max(y(i1:i2));
x_max = x_max + pos(1);

tmp = find_peaks1d(y(i1:i2), sigma, y_max/10, 1)+i1-1;

if size(tmp,1) > 1
    tmp_sort = sort([tmp y(tmp)], 2);
    p_init = double([x_max sigma y_max x(tmp_sort(2,end-1)) sigma y(tmp_sort(2,end-1)) ]);
else
    p_init = double([x_max sigma y_max x_max sigma y_max ]);
end    
%%


options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',
[param,r,J,COVB,mse] = nlinfit(x(i1:i2),  y(i1:i2) ,@two_gauss1d,p_init, options);
ci = nlparci(param,r,'Jacobian',J);
param_err = abs( ci(:,1)' - param);

param(2) = abs(param(2));
param(5) = abs(param(5));

%%
if param(1) < param(4)
    p_fast = param(4:6);
    p_fast_err = param_err(4:6);
    p_slow = param(1:3);
    p_slow_err = param_err(1:3);
else
    p_fast = param(1:3);
    p_fast_err = param_err(1:3);
    p_slow = param(4:6);
    p_slow_err = param_err(4:6);
end



%% plot results 
x_plot = x(1):(x(2)-x(1))/10:x(end);

subplot(2, 1, 1)
plot(x, y, 'g', x_plot,   gauss1d(p_slow, x_plot), 'b--', x_plot,  gauss1d(p_fast, x_plot), 'k--',  x_plot,two_gauss1d(param, x_plot), 'r--'  )
vline(x(i1), 'k--')
vline(x(i2), 'k--')
legend({'Original','slow', 'fast', 'Sum'})
set(gca, 'XLim', [x(i1)-(x2-x1) x(i2)+(x2-x1)])
set(gca, 'YLim', 1.1*[0 max(y(i1:i2))])
ylabel('Fluorescence [a.u.]')

subplot(2, 1, 2)
plot( x,   y-two_gauss1d(param, x), 'k' )
vline(x(i1), 'k--')
vline(x(i2), 'k--')
hline(0, 'k--')
set(gca, 'XLim', [x(i1)-(x2-x1) x(i2)+(x2-x1)])
%set(gca, 'YLim', 1.1*[0 max(y(i1:i2))])
xlabel('Migration Distance [pixel]')
ylabel('Residuals [a.u.]')

end

