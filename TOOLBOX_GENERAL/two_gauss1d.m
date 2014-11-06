function [ y ] = two_gauss1d( c, x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    y = c(3) .* exp( - (x-c(1)).^2 ./ (2.*c(2).^2)  ) + c(6) .* exp( - (x-c(4)).^2 ./ (2.*c(5).^2)  );

end

