function [ param ] = fit_double_band(y, I, width,  bg)
%fits a gaussian to the lane, I is intensity, y is the pixel. I(y)
%   Detailed explanation goes here

%%
options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',

p = find_peaks1d(I, width, 0.5*std(I));


i_fast = p(end);
if length(p) > 1 
    i_slow = p(end-1);
end

i_min = max(i_slow-width, 1);
i_max = min(i_fast+width, length(y) );


param_init = [y(i_slow) width I(i_slow) y(i_fast) width I(i_fast)];
    [param, chi2, residual, exitflag] = lsqcurvefit(@two_gauss1d, param_init, y(i_min:i_max), I(i_min:i_max),[],[], options);

plot(y, I, 'b', y(p), I(p), 'r.', y, gauss1d(param(1:3), y), 'g--', y, gauss1d(param(4:6), y), 'r--')
legend('data', 'maxima found', 'gaussian 1', 'gaussian 2')
pause


end

