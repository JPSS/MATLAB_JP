function [ out ] = min_max_uint16( img )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    out = img-min(img(:));
    
    out = uint16(out*(2^16-1) ./ max(out(:))    );

end

