function rgb = write2rgb(r, g, b , filepath)
%UNTITLED2 Summary of this function goes here
%   writes a 
rgb = zeros(size(r,1) , size(r,2) , 3 );

rgb(:,:,1) = r;

rgb(:,:,2) = g;

rgb(:,:,3) = b;

imwrite(uint16(rgb), filepath , 'tif');

end

