%% STARTUP
clc, clear all, close all
path0 = cd;
run('my_prefs')


%% load tform
cd(data_dir)
[map_fname map_pname]=uigetfile('*.mat','Select the mapping-file: tform_red2green.mat: ');
cd(path0)
tmp = load([map_pname map_fname] );
tform_r2g = tmp.tform;


%% select file
cd(data_dir)
[fname pname]=uigetfile('*.fits','Select the RED-movie (.fits).');
cd(path0)

path_out = pname ;

%% generate a movie file
mov_in = movie(pname, fname, 1, -1, [1]); % pname, fname, first=1, last=all, sequence=all

%% loop through movie and transform it
mov_in.initRead();
go_on = 1;
floc_out = [mov_in.pname mov_in.fname(1:end-5) '_tformed.fits'];

A_out = zeros(mov_in.sizeX, mov_in.sizeY, mov_in.mov_length, 'int16');
h = waitbar(0,'Transfroming movie... please wait');

while go_on
    [mov, frames, go_on]  = mov_in.readNext;
    
    for i=1:size(mov,3)
        A_out(:,:,frames(i)) = int16(imwarp(mov(:,:,i), tform_r2g, 'OutputView', imref2d(size(mov(:,:,i)))) ); % transfrom images
        waitbar( frames(i)/mov_in.mov_length , h, ['Transforming movie... ' num2str(frames(i)) ' of ' num2str(mov_in.mov_length) ' done']) % update waitbar

    end
    

end
fitswrite( A_out , floc_out);
close(h)

disp('Done.')
