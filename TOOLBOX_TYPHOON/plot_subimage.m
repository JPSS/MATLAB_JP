function [ output_args ] = plot_subimage(img, subimg, sigma )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if length(sigma) == 1
        sigma(2) = sigma
    end
    
    if length(subimg) < 4
        subimg = [1 1 size(img,2) size(img,1)]
        
    end
    
    colormap('Gray');
    x = reshape(img, size(img,1)*size(img,2), 1);
    imagesc(img, [mean(x)-sigma(1)*std(x) mean(x)+sigma(2)*std(x)]), axis([ subimg(1) subimg(1)+subimg(3) subimg(2) subimg(2)+subimg(4)])
    %[ subimg(2)+subimg(4) subimg(1)+subimg(3)]
end

