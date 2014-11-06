function [ z ] = generate_template( box_dim, d, s, phi )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %d=4.8; %distance

    x_m = box_dim/2;
    y_m = box_dim/2;

    dx = d*cos(phi)/2;
    dy = d*sin(phi)/2;
    
    x_a=x_m-dx; %1st peak position
    y_a=y_m-dy; %1st peak position
    x_b=x_m+dx;
    y_b=y_m+dy;    
    
   % s=1.3;  %sigma 

    z = zeros(box_dim,box_dim);


     for i=1:box_dim % i is 1pxl step
         for j=1:box_dim % j is 1pxl step
            z(i,j) = exp(-((((i-x_a)^2)+(j-y_a)^2))/(2*(s)^2)) + 2*exp(-((((i-x_b)^2)+(j-y_b)^2))/(2*(s)^2));
         end
     end

 
 
end

