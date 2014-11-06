function [ lanes ] = find_lanes_inflection( img, pos )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

area = img( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3));
x = double(pos(1):pos(1)+pos(3));


y = sum(area);

%%
kernel_sz = 21;

kernel = ones(kernel_sz,1)/kernel_sz;  % size of kernel
y_filtered = conv(y, kernel);
y_filtered = y_filtered((kernel_sz-1)/2+1:end-(kernel_sz-1)/2);
yp = abs(diff(y_filtered));
xp = x(1:end-1) + (x(2)-x(1))/2;%   x(1:length(yp));
%%
subplot(2, 1, 1)
plot(x, y, 'b', x, y_filtered, 'b--')
legend({'profile', 'filtered profile'})
set(gca, 'XLim', [pos(1) pos(1)+pos(3)])

subplot(2, 1, 2)
plot(xp, yp, 'b' ), hold on
hline(std(yp), 'r--')
xlim = [xp(1) xp(end)]; 
set(gca, 'XLim', xlim)
h = imline(gca, xlim, [std(yp) std(yp)]);
set(gca, 'XLim', [pos(1) pos(1)+pos(3)])
setColor(h,[1 0 0]);
setPositionConstraintFcn(h, @(pos)[  xlim' [pos(1,2);pos(1,2)]  ])
pos_line = wait(h);
h_min = pos_line(1,2);

%{
options.WindowStyle='normal';
prompt={'Enter threshold (default red line):'};
def={num2str(std(yp))};
threshold = inputdlg(prompt, 'Enter threshold', 1, def, options);
h_min = str2double(threshold(1));
%}

close all




%%
p = find_peaks1d(yp, kernel_sz, h_min, 1); %data, width, min_height, 1=absolute height




%%
lanes = zeros(floor(size(p,1)/2),4);



for i=1:size(lanes, 1)
   lanes(i, 2)= pos(2);
   lanes(i, 4)= pos(4);
   %lanes(i, 1) = x(p(i));
   %lanes(i, 3) = x(p(i+1))-x(p(i)); 
   
   lanes(i, 1) = pos(1)+p(2*i-1);  % xp(p(2*i-1)) wpould yield no integer
   lanes(i, 3) = p(2*i)-p(2*i-1); 

end


%{
subplot(2, 2, 2)
plot_image(img, 'Found lanes', 14);
hold on
for i=1:size(lanes,1)
    rectangle('Position', lanes(i,:), 'EdgeColor', 'r'), hold on
end

subplot(2, 2, 4)

plot(1:length(y), y, 'g', 1:length(y_filtered), y_filtered, 'r', 1:length(yp), yp, 'b'), hold on
legend({'raw', 'filtered', 'd filtered'})
for i=1:length(p)
    vline(p, 'k--')
end
%}






end

