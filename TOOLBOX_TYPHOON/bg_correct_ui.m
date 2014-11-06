function [ img_bg, bg ] = bg_correct_ui( img, img_title )
% select and subtract background from image img
%   Detailed explanation goes here

%x = reshape(img, size(img,1)*size(img,2), 1);


%colormap('Gray');
%scaling = [0 mean(x)+2*std(x)]; %[mean(x)-2*std(x) mean(x)+2*std(x)];

% imagesc(img, scaling), colorbar
%title([ img_title ' select background' ], 'FontSize' , 18)



%%
plot_image_ui(img)
title( ['Select background of ' img_title])

h = imrect;
position = wait(h);
pos = int32(getPosition(h)); % [xmin ymin width height]

bg = mean( mean(  img(pos(2):pos(2)+pos(4)  ,   pos(1):pos(1)+pos(3)   ) ));

img_bg =  img - bg ;
close all

end

 