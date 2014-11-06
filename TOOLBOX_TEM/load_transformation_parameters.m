function [ trafo ] = load_transformation_parameters( fileloc )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    data = textread(fileloc, '%s','delimiter', ';', 'emptyvalue', NaN);

    trafo = zeros(int32((size(data,1)-4)/3)+1, 13);
    for i=5:3:size(data,1)
        trafo((int32(i-4)/3)+1,:) = sscanf(data{i},'%i %i %f %f %f %f %f %f %f %f %f %f');
    end

end

