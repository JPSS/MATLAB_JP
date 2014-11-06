%%

clc
clear all
close all

[fname pname] = uigetfile('.tif', 'Select the FIRST tif image');
%%
files = dir([pname '*.tif']);
h_gauss = fspecial('gaussian', 1, 0.5) ;
%%
for i = 1:size(files,1)
    
   
   t = Tiff([pname files(i).name],'r');

   im = t.read();
      
   x = double(reshape(im, size(im,1)*size(im,2), 1));

  
   im_mean = mean(x);
   std_im = std(x);

   [n, xhist] = hist(x, 100);
   [im_max max_n] = max(n);
    
   indices = find(im > im_mean+2*std_im);
   im(indices) = xhist(max_n);
   
 %  filename =  [pname 'filtered' filesep 'filtered_' files(i).name];
 %  t = Tiff(filename,'w');
  
 
   t_out = Tiff([pname 'filtered' filesep files(i).name],'w');

   if i==1
    tagstruct.ImageLength = t.getTag('ImageLength');
    tagstruct.ImageWidth = t.getTag('ImageWidth');
    tagstruct.Photometric = t.getTag('Photometric');
    tagstruct.BitsPerSample = t.getTag('BitsPerSample');
    tagstruct.SamplesPerPixel = t.getTag('SamplesPerPixel');
    tagstruct.RowsPerStrip = t.getTag('RowsPerStrip');
    tagstruct.PlanarConfiguration = t.getTag('PlanarConfiguration');
    tagstruct.Software = 'MATLAB';
   end
    
    t_out.setTag(tagstruct)
   
   t_out.write(im)
   
  %imwrite(im, 'tif')
   
   if mod(i,50) == 0
       display(num2str(i));
   end
   
   
   
end

display('done')
%%don

