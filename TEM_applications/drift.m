%%
run('my_prefs')

%%

files = dir(cd)

%%
images = cell(5,1);

for i=3:7
    images{i-2} = imread(files(i).name);
end
%%

img = zeros(1024,1024,5);

for i=1:5
    img(:,:,i) = images{i}(1:1024, 1:1024);
end
    
%%
img = uint16(img);

%%
cor11 = normxcorr2(img(:,:,1), img(:,:,1)); 
cor12 = normxcorr2(img(:,:,1), img(:,:,2)); 
cor23 = normxcorr2(img(:,:,2), img(:,:,3));
cor34 = normxcorr2(img(:,:,3), img(:,:,4)); 
%%
close all

subplot(1,2,1)
imagesc(cor11), axis image, colorbar

subplot(1,2,2)

%%
close all
imagesc(cor12), axis image, colorbar

%%
%cor23 = normxcorr2(img(:,:,2), img(:,:,3)); 
close all
imagesc(cor23), axis image, colorbar

%%
%cor34 = normxcorr2(img(:,:,3), img(:,:,4)); 
close all
imagesc(cor34), axis image, colorbar





