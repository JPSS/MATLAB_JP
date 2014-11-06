%% startup
clear all, close all%, clc
run('my_prefs')
path0=cd;

% load dm4 file
%cd(data_dir)
cd('/Users/jonasfunke/Documents/FRET_STAGE/TEM images/2013-11-05_cryo-data')
[fname pname] = uigetfile('.dm4', 'Select cryo-em average image');
cd(path0)

%
[m , sx, units] = ReadDMFile([pname fname]);
display('Image(s) loaded')
% write image to 
%m = m-min(m(:));
%imwrite( uint16(m.*(2^16-1)./max(m(:))), [pname fname(1:end-4) '.tif'])
%%
path_out = pname; %[pname fname(1:end-4) '_uint16_scaled'];
%mkdir(path_out)


left = min(m(:));
right = max(m(:)-left);
imwrite( uint16(   (m-left).*(2^16-1)./right), [path_out fname(1:end-4) '.tif'])

m_inv = -m;
left = min(m_inv(:));
right = max(m_inv(:)-left);
imwrite( uint16(   (m_inv-left).*(2^16-1)./right), [path_out fname(1:end-4) '_inverted.tif'])





%{

%%
for i=1:size(m,3)
    t = Tiff([path_out filesep fname(1:end-4) '_' sprintf('%.2i', i) '.tif'],'w');
    t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
    t.setTag('BitsPerSample',32);
    t.setTag('SampleFormat',Tiff.SampleFormat.IEEEFP);
    t.setTag('ImageLength',size(m,1));
    t.setTag('ImageWidth',size(m,2));
    t.setTag('SamplesPerPixel',1);
    t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    t.write(  m(:,:,i)  );
    t.close();
end
%%
path_out2 = [pname fname(1:end-4) '_32bit_scaled'];
mkdir(path_out2)

left = min(m(:));
right = max(m(:)-left);

for i=1:size(m,3)
    t = Tiff([path_out2 filesep fname(1:end-4) '_' sprintf('%.2i', i) '.tif'],'w');
    t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
  %  t.setTag('Compression',Tiff.Compression.None);
    t.setTag('BitsPerSample',32);

    t.setTag('SamplesPerPixel',1);
  %  t.setTag('SampleFormat',Tiff.SampleFormat.IEEEFP);
   % t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);

    t.setTag('ImageLength',size(m,1));
    t.setTag('ImageWidth',size(m,2));
  %  t.setTag('TileLength',32);
  %  t.setTag('TileWidth',32);
    t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    
    t.write(  uint32((m(:,:,i)-left).*(2^32-1)./right)  );
    %t.write(  m(:,:,i)  );
    
    t.close();
end


%%


if size(m,3)>1
    path_out_subavg = [path_out filesep 'subavg'];
    mkdir(path_out_subavg)
    subavg = single(zeros(size(m,1), size(m,2), 7));
    for i=1:7
        tmp = zeros(size(m,1), size(m,2));
        for j=0:3
            tmp = tmp + m(:,:,i+j);
        end

        subavg(:,:,i) = single(tmp ./4);
        %imwrite( uint16(  subavg(:,:,i)    ), [path_out filesep fname(1:end-4) '_' sprintf('%.2i', i) '_' num2str(i) '-' num2str(i+j) '.tif']) % no normalization

        t = Tiff([path_out_subavg filesep fname(1:end-4) '_' sprintf('%.2i', i) '_subavg_' num2str(i) '-' num2str(i+j) '.tif'],'w');
        t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
        t.setTag('BitsPerSample',32);
        t.setTag('SampleFormat',Tiff.SampleFormat.IEEEFP);
        t.setTag('ImageLength',size(subavg,1));
        t.setTag('ImageWidth',size(subavg,2));
        t.setTag('SamplesPerPixel',1);
        t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
        t.write(  subavg(:,:,i)  );
        t.close();
    end
    

    t = Tiff([path_out filesep fname(1:end-4) '_avg.tif'],'w');
    t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
    t.setTag('BitsPerSample',32);
    t.setTag('SampleFormat',Tiff.SampleFormat.IEEEFP);
    t.setTag('ImageLength',size(subavg,1));
    t.setTag('ImageWidth',size(subavg,2));
    t.setTag('SamplesPerPixel',1);
    t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    t.write(  single(sum(double(m),3)./size(m,3))  );
    t.close();

end

%}
display(['Done writing: ' fname])


