function [lane, y, pos] = integrate_lane( img, img_title )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



    cur_fig = plot_image_ui(img)  ;
    title(img_title)
    h = imrect;
    position = wait(h);
    pos = int32(getPosition(h)); % [xmin ymin width height]
    y = transpose(pos(2):pos(2)+pos(4));
    lane = transpose(sum(transpose(  img(y  ,   pos(1):pos(1)+pos(3)   )  )) );

    close(cur_fig)
    
   

end

