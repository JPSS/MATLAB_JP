function rgb = write2color(im,  filepath, color)
%writes an image im to a rgb tif
%   writes a 
rgb = zeros(size(im,1) , size(im,2) , 3 );

if strcmpi(color,'red')
    rgb(:,:,1) = im;
end

if strcmpi(color,'green')
    rgb(:,:,2) = im;
end

if strcmpi(color,'blue')
    rgb(:,:,3) = im;
end

imwrite(uint16(rgb), filepath , 'tif');

end

