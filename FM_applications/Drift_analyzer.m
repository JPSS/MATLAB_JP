%% STARTUP
clc, clear all, close all
path0 = cd;
run('my_prefs')

%% load movie
cd(data_dir)
[fname pname]=uigetfile('*.fits','Select a movie (.fits).');
cd(path0)

%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', 'Sequence:' };
input_default = {'2', '-1', '1'};
tmp = inputdlg(input, 'Parameters', 1, input_default);


first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file
%determine sequences 
sequence = zeros(1, size(tmp{3},2));
for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence(1,i) =1;
    end
end


%% init movie class
mov = movie(pname, fname, first, last, sequence); % pname, fname, first=1, last=all, sequence=all

%% compute xcorrelation

xcorr = zeros(mov.sizeX, mov.sizeY, mov.mov_length);

h = waitbar(0,'Calculating xcorrelation... please wait');

corr_trace = zeros(length(mov.frames),3);

go_on = 1;
mov.initRead();
N = 1;
while go_on
    [cur_mov, frames, go_on]  = mov.readNext;
    
    if N==1
        first_img = cur_mov(:,:,1);
    end
    
    for i=1:size(cur_mov,3)

        tmp = normxcorr2(cur_mov(:,:,i), first_img);
        
        [v,ind]=max(tmp(:));
        [corr_trace(N,2), corr_trace(N,3)] = ind2sub(size(tmp),ind);
        corr_trace(N,1)  = v;
        N = N+1;

        waitbar( frames(i)/mov.frames(end) , h, ['Calculating xcorrelation... ' num2str(frames(i)) ' of ' num2str(mov.frames(end)) ' done']) % update waitbar
        
    end
    

end
close(h)

%% plot dx and dx against frames
xy_corr = [corr_trace(:,2)-corr_trace(1,2) corr_trace(:,3)-corr_trace(1,3) ] ;
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(2, 1, 1)
plot(mov.frames, corr_trace(:,1), 'k', 'Linewidth', 1)
xlabel('Frame'), ylabel('Correlation Coefficient')
set(gca, 'YLim', [0 1])

subplot(2, 1, 2)
plot(mov.frames, xy_corr(:,1), 'g', mov.frames, xy_corr(:,2), 'b', 'Linewidth', 1)
legend({'x_i-x_1', 'y_i-y_1'})
xlabel('Frame'), ylabel('Dirift [pixel]')
ylim = [min(xy_corr(:))-1 max(xy_corr(:))+1];
set(gca, 'YLim', ylim)
print(cur_fig, '-dtiff', '-r300', [pname filesep fname(1:end-5) '_drift.tif'])

%%
disp('Done.')