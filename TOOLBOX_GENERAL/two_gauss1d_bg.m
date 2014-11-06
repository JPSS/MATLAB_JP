function [ y ] = two_gauss1d_bg( c, x )
% fits the sum of two gaussian an background
%   Detailed explanation goes here

    y = c(3) .* exp( - (x-c(1)).^2 ./ (2.*c(2).^2)  ) + c(6) .* exp( - (x-c(4)).^2 ./ (2.*c(5).^2)  ) + c(7);

end