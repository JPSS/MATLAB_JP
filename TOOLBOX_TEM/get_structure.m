function [ structure, tmp ] = get_structure( file_loc, name, length )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
structure = cell(1,4);
tmp = dlmread(file_loc, '\t',1); % skip first line
disp('skip first line')
    % tmp = dlmread(file_loc, '\t',1);


%% check for double structre and so on
dangle = zeros(size(tmp,1)/2, 1);
slice = zeros(size(tmp,1)/2, 1);
for i=1:2:size(tmp,1)
        if tmp(i,7) == tmp(i+1,7)   
            dangle((i+1)/2) = abs(tmp(i,6) - tmp(i+1,6)) ;
            if dangle((i+1)/2) > 180
                 dangle((i+1)/2) = 360 - dangle((i+1)/2);
            end
            slice( (i+1)/2 ) = tmp(i,7);
        else
            display(['WARNING: Mismatch, ' num2str(tmp(i,7)) ' not equal to ' num2str(tmp(i+1,7)) '. Check data table!']);
            i
        end
end
structure{1,1} = dangle ;
structure{1, 2} = slice;
structure{1, 3} = length; 
structure{1, 4} = name;

%%
end

