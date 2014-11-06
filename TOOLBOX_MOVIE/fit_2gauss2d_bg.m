function [ result, err, ci, area, residuals ] = fit_2gauss2d_bg( x0, y0, x1, y1, sigma, w_fit, img ) % (x0, y0, s, w, image)
%UNTITLED2 Summary of this function goes here
%   c_init = [x_0 y_0 s_x s_y A-bg bg];
    
    options = optimset('Algorithm','levenberg-marquardt','display','final', 'MaxFunEvals',50000,'TolFun',1e-9,'MaxIter',200, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',

    xm = round((x0+x1)/2)
    ym = round((y0+y1)/2)
    x0 = round(x0);
    x1 = round(x1);
    y0 = round(y0);
    y1 = round(y1);
    
    area = [  max(xm-w_fit,1) , min(xm+w_fit,size(img,1)) , max(ym-w_fit,1) , min(ym+w_fit,size(img,2)) ];
    [X,Y] = meshgrid(area(1):area(2) , area(3):area(4) );
    Z = img(area(3):area(4), area(1):area(2));
    xyz_data=zeros(size(X,1)*size(X,2),3);
    xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
    xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix X
    xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix X
    
    
    % initial parameters
    bg = mean(xyz_data(:,3));
    %c_init = [x0 y0 sigma img(y0,x0)-bg x1 y1 sigma img(y1,x1)-bg bg]
    c_init = [x0 y0 sigma img(y0,x0)-bg x1 y1 img(y1,x1)-bg bg];
    
    same_width_gauss = @(c,x) two_gauss2d_bg([c(1:6) c(3) c(7:8)], x);
    
    
    % fit the parameters
    [coef,residuals,J,COVB,mse, ErrorModelInfo] = nlinfit(  xyz_data(:,1:2)  ,  xyz_data(:,3), same_width_gauss, c_init, options);
    
    
    coef(3) = abs(coef(3));
    %coef(7) = abs(coef(7));

    % calculate errors
    ci = nlparci(coef,residuals,'Jacobian',J);
    
    % caculate c_err
    err = coef-ci(:,1)'+ci(:,2)'-coef;
    
    
    err = [err(1:6) err(3) err(7:8)];
    ci = [ci(1:6) ci(3) ci(7:8)];

    result = [coef(1:6) coef(3) coef(7:8)];
    
    disp(ErrorModelInfo)
    
end

