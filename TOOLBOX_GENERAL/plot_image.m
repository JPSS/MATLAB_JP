function [ ] = plot_image( img, frac )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dc = max(img(:))-min(img(:));
    if length(frac)>1
        clim = [median(img(:))-frac(1)*dc median(img(:))+frac(2)*dc];
    else
        clim = [median(img(:))-frac*dc median(img(:))+frac*dc];
    end
    imagesc(img, clim), axis image, colormap gray
end

