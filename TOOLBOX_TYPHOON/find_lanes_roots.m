function [ lanes ] = find_lanes_roots( img, pos )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

area = img( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3));
x = double(pos(1):pos(1)+pos(3));



y = sum(area);

%%
close all
subplot(2, 1, 1)
imagesc(img), axis image, colormap gray
set(gca, 'XLim', [pos(1) pos(1)+pos(3)])
set(gca, 'YLim', [pos(2) pos(2)+pos(4)])

shift = mean(y)-std(y);
subplot(2, 1, 2)
hleg(1) = plot(x, y, 'b');
hold on
hleg(2) = hline(mean(y), 'k--');
hleg(3) = hline(min(y), 'y--');
hleg(4) = hline(max(y), 'g--');
y_init = (max(y)-min(y))/2;
hleg(5) = hline(y_init, 'r--');

legend(hleg, {'y', 'mean(y)', 'min(y)', 'max(y)', '(max(y)-min(y))/2 (default)'})
set(gca, 'XLim', [pos(2) pos(2)+pos(4)])

plot(x, y, 'b')
xlim = [x(1) x(end)]; 
set(gca, 'XLim', xlim)
h = imline(gca, xlim, [y_init y_init]);
setColor(h,[1 0 0]);
setPositionConstraintFcn(h, @(pos)[  xlim' [pos(1,2);pos(1,2)]  ])
pos_line = wait(h);
shift = pos_line(1,2);
%{
options.WindowStyle='normal';
prompt={'Enter shift (default red line):'};
def={num2str(mean(y)-std(y))};
threshold = inputdlg(prompt, 'Enter shift', 1, def, options);
shift = str2double(threshold(1));
close all
%}

close all
%%
y_shift = y - shift;

%find roots of y_shift
r = [];
for i=1:length(y_shift)-1
    if y_shift(i)*y_shift(i+1) < 0 %root
       r = [r ; i];
    end
end





lanes = zeros(length(r)/2,4);



for i=1:size(lanes, 1)
   lanes(i, 2)= pos(2);
   lanes(i, 4)= pos(4);
   %lanes(i, 1) = x(p(i));
   %lanes(i, 3) = x(p(i+1))-x(p(i)); 
   
   lanes(i, 1) = x(r(i*2-1));
   lanes(i, 3) = x(r(i*2))-x(r(i*2-1)); 

end


%{
subplot(2, 2, 1)
plot_image(img, 'Found lanes', 14);
hold on
for i=1:size(lanes,1)
    rectangle('Position', lanes(i,:), 'EdgeColor', 'r'), hold on
end

subplot(2, 2, 3)
plot(x, y, 'b', x, y_shift, 'g', x(r), y(r), 'r.')
legend('y', 'y-y_mean-y_std', 'roots')
%}




end

