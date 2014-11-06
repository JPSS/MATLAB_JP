function [ lanes , area] = find_lanes( img )
% autmatically find lanes based on inflection-point-algorithm or simple
% threshold (roots)
% user can choose which result to use
more = 1;
lanes = [];
area = [];
while more
    close all


    %% plot image
    plot_image_ui(img)
    title('Select area of lanes')
    h = imrect;
    position = wait(h);
    pos = int32(getPosition(h)); % [xmin ymin width height]

    area = [area ; pos];
    %% 1st alogrithm, uses the roots
    lanes_roots = find_lanes_roots(img, pos);

    %% 2nd algorithm, uses inflection points
    lanes_inflection = find_lanes_inflection(img, pos);

    %% show results and ask which result to use
    close all
    scrsz = get(0,'ScreenSize');
    figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)/2]);

    subplot(1, 2, 1)
    imagesc(img), axis image, colormap gray
    title('Root algorithm');
    hold on
    for i=1:size(lanes_roots,1)
        rectangle('Position', lanes_roots(i,:), 'EdgeColor', 'r'), hold on
    end

    subplot(1, 2, 2 )
    imagesc(img), axis image, colormap gray
    title('Infelction-point algorithm');
    hold on
    for i=1:size(lanes_inflection,1)
        rectangle('Position', lanes_inflection(i,:), 'EdgeColor', 'r'), hold on
    end

    button = questdlg('Use left or right lanes?','lane-finding algorithm' ,'Left','Right', 'Left');
    if strcmp(button,'Left')
        lanes = [lanes ; lanes_roots];
    else
        lanes = [lanes ; lanes_inflection];
    end

    button = questdlg('Do you want to add more lanes?','more lanes' ,'Add Lanes','Enough', 'Enough');
    if strcmp(button,'Add Lanes')
        more = 1;
    else
        more = 0;
    end

    close all
end

display(['Found '  num2str(size(lanes,1)) ' lanes.'])

end







%% old stuff


%{

%%
nlanes = 20;
w = pos(3) /(4*nlanes)
wsz = double(max(1,w/2))

y = filter(ones(1,wsz) ./ wsz,1,sum(a));
h_min = std(y);


p = find_peaks1d(y, w, h_min);

wsz = double(max(1,w/4))
dy = diff(y);

dy_filtered = filter( ones(1,wsz) ./ wsz,1, dy);
ddy = diff(dy);
ddy_f = diff(dy_filtered);


%{
subplot(2, 1, 1)

plot(dy, 'b'), hold on

plot(dy_filtered, 'r')

subplot(2, 1, 2)

plot(ddy, 'b'), hold on
plot(ddy_f, 'r')
%}


%%
p_ddy_f = find_peaks1d(-ddy_f, w/2, 0.5*std(-ddy_f));

%% plot everything
%lanes = zeros(length(p)-1,4);

lanes = zeros(length(p_ddy_f)/2,4);



for i=1:size(lanes, 1)
   lanes(i, 2)= pos(2);
   lanes(i, 4)= pos(4);
   %lanes(i, 1) = x(p(i));
   %lanes(i, 3) = x(p(i+1))-x(p(i)); 
   
   lanes(i, 1) = x(p_ddy_f(i*2-1));
   lanes(i, 3) = x(p_ddy_f(i*2))-x(p_ddy_f(i*2-1)); 

end




subplot(3, 1, 1)
plot_image(img, 'Found lanes', 14);
hold on
for i=1:size(lanes,1)
    rectangle('Position', lanes(i,:), 'EdgeColor', 'r'), hold on
end

subplot(3, 1, 2)
plot(x, -y, 'b', x(p) , -y(p), 'r.');

subplot(3,1, 3)

plot(x(1:end-2), -ddy_f, 'b', x(p_ddy_f), -ddy_f(p_ddy_f), 'r.')



%}
