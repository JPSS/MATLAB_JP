clear all, close all, clc
path0 = cd;
%% Load ml2d.doc and get parameters
[fname pname] = uigetfile('*.doc','Select ml2d.doc file ');
trafo = load_transformation_parameters([pname fname]);

%% select folder wwith all tif images
[fname_first pname_first] = uigetfile('*.tif','Select first file ');

cd(pname_first)
[fname_last pname_last] = uigetfile('*.tif','Select first file ');
cd(path0)
%%
start = str2double(fname_first(end-6:end-4))
stop = str2double(fname_last(end-6:end-4))
n_img = stop-start+1

%% loop though images and rotate them
close all
for i=1:100%size(trafo,1)
    
    if trafo(i,8)==8
        trafo(i,:)
        img = double(imread([pname_last filesep fname_first(1:end-6) sprintf('%.2i', trafo(i,1)-1) '.tif']));

        scale = [mean(img(:))-2*std(img(:)) mean(img(:))+2*std(img(:))];

        subplot(1, 3, 1)
        imagesc(img, scale), colorbar, colormap gray, axis image
        title(['Class ' num2str(trafo(i,8))])

        if trafo(i,9)==1
            img_trafo = flipdim(img,1);
        else
            img_trafo = img;
        end
        subplot(1,3,2)
        alpha =  trafo(i,5) ;
        imagesc(img_trafo, scale), colorbar, colormap gray, axis image
        
                subplot(1,3,3)
        alpha =  trafo(i,5) ;
        imagesc(imrotate(img_trafo, -alpha, 'nearest', 'crop'), scale), colorbar, colormap gray, axis image
        
        pause
    end
    
    
end



%%


