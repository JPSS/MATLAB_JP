function [ leak_dir] = calculate_corrections(dd_bg, da_bg, aa_bg, file_out )
%determines leakage and direct excitation correction from dd_img, da_img
%and aa_img
%   Detailed explanation goes here

%%



%% fret lanes
display('Select donor-only trace!')
donly = cell(1, 7);
[l1, y , pos_d] = integrate_lane(dd_bg, 'D->D image: Select D-ONLY band');
donly{1, 1} = +1.* l1;
donly{1, 2} =  +1.*   transpose( sum(transpose(  da_bg(y  ,   pos_d(1):pos_d(1)+pos_d(3)   )  )) );
donly{1, 3} =  +1.*  transpose( sum(transpose(  aa_bg(y  ,   pos_d(1):pos_d(1)+pos_d(3)   )  )) );
donly{1, 4} = y;
donly{1,5} = pos_d;
donly{1,6} = 'D-only';

%%

display('Select acceptor-only trace!')
aonly = cell(1, 7);
[l1, y , pos] = integrate_lane(aa_bg, 'a-only: A -> A');
aonly{1, 3} = +1.* l1;
aonly{1, 2} = +1.*    transpose( sum(transpose(  da_bg(y  ,   pos(1):pos(1)+pos(3)   )  )) );
aonly{1, 1} = +1.*    transpose( sum(transpose(  dd_bg(y  ,   pos(1):pos(1)+pos(3)   )  )) );
aonly{1, 4} = y;
aonly{1,5} = pos;
aonly{1,6} = 'A-only';

%% show selction
%plot_image_ui(da_bg), hold on 
%rectangle('Position',donly{1, 5}, 'EdgeColor', [0 1 0])
%rectangle('Position',aonly{1, 5}, 'EdgeColor', [1 0 0])


%%
%{
cur_fig = figure('Visible','on', 'PaperPositionMode', 'auto'); % figure('Visible','off');%left bottom width height
x_raw = abs(donly{1}); y_raw = abs(donly{2});
x = []; y = [];
index = [];
max_x = max(x_raw);
for i=1:length(y_raw)
   if x_raw(i) >= 0.4*max_x
      x = [x ; x_raw(i) ];
      y = [y ; y_raw(i) ];
   end
end
p_leak = polyfit(x, y, 1);

subplot(1, 2 , 1)
plot(x_raw, y_raw, 'b.', x, y, 'g.', x , p_leak(1)*x+p_leak(2), 'r', 'MarkerSize', 10 )
legend('data', 'data > 0.4*max', 'fit')
xlabel('I_{D->D}')
ylabel('I_{D->A}')
title('Leakage correction')


x_raw = abs(aonly{3}); y_raw =abs(aonly{2});
x = []; y = [];
index = [];
max_x = max(x_raw);
for i=1:length(y_raw)
   if x_raw(i) >= 0.4*max_x
      x = [x ; x_raw(i) ];
      y = [y ; y_raw(i) ];
   end
end
p_dir = polyfit(x, y, 1);

subplot(1, 2, 2)
plot(x_raw, y_raw, 'b.', x, y, 'g.', x , p_dir(1)*x+p_dir(2), 'r', 'MarkerSize', 10)
legend('data', 'data > 0.4*max', 'fit')
xlabel('I_{A->A}')
ylabel('I_{D->A}')
title('Direct-excitation correction')



print(cur_fig, '-djpeg' , '-r300', [file_out(1:end-4) '_leak_dir_fit.jpg']); %save figure




%}


%%
pos = donly{5};

lane_dd = dd_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
lane_da = da_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));

[cc, shift, lane_da_shift] = xcorr2_bounded(lane_dd, lane_da, 5, 0); % find best overlay of images

x_raw = reshape(lane_dd, size(lane_dd,1)*size(lane_dd,2), 1);
y_raw = reshape(lane_da_shift, size(lane_da_shift,1)*size(lane_da_shift,2), 1);

mean_x = mean(x_raw);
mean_y = mean(y_raw);

i_found = find( x_raw > mean_x & y_raw > mean_y); % only use intensities wich are greater than mean (dont fit to background)
x = x_raw(i_found);
y = y_raw(i_found);

p_leak = polyfit(x, y, 1);

scrsz = get(0,'ScreenSize');
fig_dim =[2*13.7 2*8];%5.23];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(2, 2, 1)
xlim = [min(x_raw) max(x_raw)];
plot(x_raw, y_raw, 'b.', x, y, 'g.', xlim , p_leak(1)*xlim+p_leak(2), 'r--', 'MarkerSize', 1 )
legend({'Discarded', 'Used to fit', ['Fit, dir=' num2str(round(p_leak(1)*100)) '%']}, 'FontSize', 10, 'Location', 'SouthEast')
xlabel('D->D')
ylabel('D->A')
title('Leakage correction')
set(gca, 'XLim', [xlim(1) xlim(2)*1.2])

subplot(2, 2, 2)
plot(donly{4}, donly{1}, 'g', donly{4}, donly{2}, 'b', donly{4}, donly{2}-donly{1}*p_leak(1,1), 'b--' )
legend({'D->D', 'D->A', 'D->A corrected'}, 'FontSize', 10)
%title('D-Only')
xlabel('Migration distance [pixel]')
ylabel('Intensity')
set(gca, 'XLim', [donly{4}(1) donly{4}(end)])




%
pos = aonly{5};

lane_aa = aa_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
lane_da = da_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));

[cc, shift, lane_da_shift] = xcorr2_bounded(lane_aa, lane_da, 5, 0); % find best overlay of images

x_raw = reshape(lane_aa, size(lane_aa,1)*size(lane_aa,2), 1);
y_raw = reshape(lane_da_shift, size(lane_da_shift,1)*size(lane_da_shift,2), 1);

mean_x = mean(x_raw);
mean_y = mean(y_raw);

i_found = find( x_raw > mean_x & y_raw > mean_y); % only use intensities wich are greater than mean (dont fit to background)
x = x_raw(i_found);
y = y_raw(i_found);

p_dir = polyfit(x, y, 1);

subplot(2,2,3)
xlim = [min(x_raw) max(x_raw)];
plot(x_raw, y_raw, 'b.', x, y, 'g.', xlim , p_dir(1)*xlim+p_dir(2), 'r--', 'MarkerSize', 1 )
legend({'Discarded', 'Used to fit', ['Fit, dir=' num2str(round(p_dir(1)*100)) '%']}, 'FontSize', 10, 'Location', 'SouthEast')
xlabel('A -> A')
ylabel('D -> A')
%title('Leakage correction')
set(gca, 'XLim', [min(x_raw) max(x_raw)*1.2])
%print(cur_fig, '-depsc2','-loose' , [file_out(1:end-4) '_dir_1']); %save figure

subplot(2,2,4)
plot(aonly{4}, aonly{3}, 'r', aonly{4}, aonly{2}, 'b', aonly{4}, aonly{2}-aonly{3}*p_dir(1), 'k--' )
%plot(donly{1}, donly{4}, 'g', donly{2}, donly{4}, 'b',  donly{2}-donly{1}*p_leak(1,1), donly{4}, 'b--' )

legend({'A->A', 'D->A', 'D->A corrected'}, 'FontSize', 10)
%title('D-Only')
xlabel('Gel-Koordinate [pixel]')
ylabel('Intesit\"at')
set(gca, 'XLim', [aonly{4}(1) aonly{4}(end)])
print(cur_fig, '-depsc2','-loose' , [file_out(1:end-4) '_corrections']); %save figure


%display('Press some key to continue...')
questdlg('Go on','Halt','Go on','Go on');

close all


%% compine to matrix
leak_dir = [p_leak ; p_dir]; 
disp(['Leakage: ' num2str(round(p_leak(1)*100)) '% Direct Ex.: ' num2str(round(p_dir(1)*100)) '%' ])

%% save data
save(file_out, 'leak_dir' ,'-ascii')

close all











end

