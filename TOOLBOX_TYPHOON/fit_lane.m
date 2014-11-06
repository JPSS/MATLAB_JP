function [ param param_err] = fit_lane(y, I, width, sigma_left, sigma_right, bg)
%fits a gaussian to the lane, I is intensity, y is the pixel. I(y)
%   Detailed explanation goes here
options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',

%y=transpose(1:max(size(I)));

%deternime initial paramteters
[height meanI] = max(I);   %maximum of lane, heigth of peak
i_min = max(1, meanI-width);    
i_max = min(meanI+width, length(y));

%fit to all data
if bg == 0  %no bg 
    param_init = [y(meanI) width height];
    
    %keyboard
    [param1, chi2, residual, exitflag] = lsqcurvefit(@gauss1d, param_init, y(i_min:i_max), I(i_min:i_max),[],[], options);
 
    %use std to reduce fitting fit again
    i_min = max(round(meanI-sigma_left*abs(param1(2))), 1 ) ;
    i_max = min(round(meanI+sigma_right*abs(param1(2))), length(y));

    param_init = param1(1:3);% param_init(4)=0;
   % [param, chi2, residual, exitflag] = lsqcurvefit(@gauss1d, param_init, y(i_min:i_max), I(i_min:i_max),[],[], options);
    
    [param,r,J,COVB,mse] = nlinfit(y(i_min:i_max),  I(i_min:i_max) ,@gauss1d,param_init, options);
    ci = nlparci(param,r,'Jacobian',J);
    param_err = abs( ci(:,1)' - param);
    
    param = [param 0 ];
        param_err = [param_err 0 ];

else %with bg
    param_init = [y(meanI) width height-I(end) I(end)];
    [param1, chi2, residual, exitflag] = lsqcurvefit(@gauss1d_bg, param_init, y(i_min:i_max), I(i_min:i_max),[],[], options);
 
    %use std to reduce fitting fit again
    i_min = max(round(meanI-sigma_left*abs(param1(2))), 1 ) ;
    i_max = min(round(meanI+sigma_right*abs(param1(2))), length(y));

    param_init = param1; 
    [param, chi2, residual, exitflag] = lsqcurvefit(@gauss1d_bg, param_init, y(i_min:i_max), I(i_min:i_max),[],[], options);
    
    
end
param(2) = abs(param(2));



%plot(y, I, 'b', y, gauss1d_bg(param1, y), 'r',y, gauss1d_bg(param, y), 'g--')
%title(['Fit from ' num2str(y(i_min)) ' to ' num2str(y(i_max))])
% legend('Data', '1st fit', '2nd fit')


end

