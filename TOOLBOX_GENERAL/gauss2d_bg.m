function [ z ] = gauss2d_bg(parameter, x) %x = (x_coor, y_coor)
%gauss2d_bg computes a vector z of each high such as x, y, z triples
    x_0 = parameter(1);
    y_0 = parameter(2);
    s_x = parameter(3);
    s_y = parameter(4);
    A = parameter(5);
    bg = parameter(6);
    z = A * exp(- ( (x(:,1)-x_0).^2./(2*s_x^2) + (x(:,2)-y_0).^2./(2*s_y^2) ) ) + bg;

end

