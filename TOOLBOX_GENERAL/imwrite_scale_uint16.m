function [ ] = imwrite_scale_uint16( img, location )
% scales image to full uni16 range and writes it
    img = img - min(img(:));
    imwrite(uint16(img .* (2^16-1) ./ max(img(:))), location)
end

