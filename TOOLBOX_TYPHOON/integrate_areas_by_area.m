function [ I, areas ] = integrate_areas_by_area( img, n_areas, area )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
   close all
    areas = zeros(n_areas, 4);
    I = zeros(n_areas, max(size(img)));
    plot_image_ui(img{1})
    
    for i=1:n_areas
        for j=1:i-1
            rectangle('Position', areas(j,:), 'EdgeColor', 'r'); %plot all integrated areas
        end

        h = imrect(gca, double(area));
        setResizable(h, 0)
        wait(h);
        pos = int32(getPosition(h)); % [xmin ymin width height]
        delete(h)
        areas(i,:) = pos;
        for j= 1:max(size(img))
            I(i,j) = sum(sum((img{j}(pos(2):pos(2)+pos(4)  , pos(1):pos(1)+pos(3))))); %integrate
        end
    end
    
    close all

end

