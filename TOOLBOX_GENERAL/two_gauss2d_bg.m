function [ z ] = two_gauss2d_bg(parameter, x) %x = (x_coor, y_coor)
%gauss2d_bg computes a vector z of each high such as x, y, z triples
    x_0 = parameter(1);
    y_0 = parameter(2);
    s_0 = parameter(3);
    A = parameter(4);
    
    x_1 = parameter(5);
    y_1 = parameter(6);
    s_1 = parameter(7);
    B = parameter(8);    
    bg = parameter(6);
    z = A * exp(- ( (x(:,1)-x_0).^2./(2*s_0^2) + (x(:,2)-y_0).^2./(2*s_0^2) ) ) + B * exp(- ( (x(:,1)-x_1).^2./(2*s_1^2) + (x(:,2)-y_1).^2./(2*s_1^2) ) ) + bg;

end

