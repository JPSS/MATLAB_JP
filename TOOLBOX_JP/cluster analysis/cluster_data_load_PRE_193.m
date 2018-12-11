function [ results, mean_times ] = cluster_data_load_PRE_193(nr_structures,nr_time_steps, varargin)
%% analyzes gillespie simulation data
%   nr_structures is number of structures in file
%   nr_times_steps is number of steps printed
%   structureNames is cell of supplied structure names

%   filename{N} are filenames of files
%   pathname is pathnames of files

%% parse input variables
    p = inputParser;
    % required parameter
    addRequired(p,'nr_structures');
    addRequired(p,'nr_times_steps');
    
%% load data
%opens file, copies data into results and currentData, analyzes,
%then opens next file
%current_data is array of averaged traces from file

[filename, pathname]=uigetfile('*','select cluster data','MultiSelect','on');

nr_parameter_lines=28;                                               %number of lines after structure names before data starts
parameters_text=cell(1+nr_structures+nr_parameter_lines,1);           %saves parameters text, each cell one line
nr_parameters=10;                                                   %number of parameters in filename

nr_files=size(filename,2);                                          %number of files loaded

if ~(iscell(filename))                                              %check if a single file was loaded, if yes, adjust data type of filename
    nr_files=1;
    filename={filename};
end

%results=cell(numFiles+1,1+nr_parameters+3*nr_structures);                         %analysis results of data
%averagedData=zeros(nr_times_steps,nr_structures,numFiles);

for file_nr=1:nr_files
    file_nr
    current_file=fopen([pathname filename{file_nr}],'r');                  %load current file
    filename{file_nr}                                                     %print filename
    for j=1:1+nr_structures+nr_parameter_lines
        parameters_text{j}=fgetl(current_file);                           %load parameter text
    end
   
    current_data=zeros(nr_time_steps+1,24,nr_structures);
    
    %temporary storage current mean folding times
    current_mean_times = zeros(nr_structures, 10);
    
%     current_line = fgetl(current_file)
%     current_line = fgetl(current_file)
%     current_line = fgetl(current_file)

    
    current_line = fgetl(current_file);
    trace_nr = zeros(nr_structures,1);

    %load simulation data line by line
    while and(current_line ~= -1, length(current_line) ~= 1)
        
        %set current_time to max possible value = max simulation duration
        current_times = nr_time_steps * ones(nr_structures, 10);
        
        current_structure = sscanf(current_line,'%*s %d',[1,2]);    %determine which structure trace belongs to
        current_structure = current_structure(1)+1;
        for step_nr=1:1+nr_time_steps
            current_line=fgetl(current_file);
            current_line_data = sscanf(current_line,'%f',[1,24]);
            current_data(step_nr,:,current_structure)=current_data(step_nr,:,current_structure) + current_line_data;
            
            %set current folding time to time at which good crossover ratio is > cutoff fraction
            for cutoff = 1:10
                if current_line_data(14) > cutoff/10 && current_times(current_structure, cutoff) == nr_time_steps
                    current_times(current_structure, cutoff) = step_nr;
                end
            end
        end
        %add current folding times to current mean times
        current_mean_times(current_structure, :) = current_mean_times(current_structure, :) + current_times(current_structure, :);
        
        current_line = fgetl(current_file);
        trace_nr(current_structure) = trace_nr(current_structure) + 1;
    end   

    fclose(current_file);
    for current_structure=1:nr_structures
        mean_times(((file_nr - 1) * nr_structures ) + current_structure, :) = current_mean_times(current_structure, :) / trace_nr(current_structure);
        results{((file_nr - 1) * nr_structures ) + current_structure ,1} = current_data(:,:,current_structure)/trace_nr(current_structure);
    end

    %load parameter valuestart and end positions from file name
    [parameterStarts, parameterEnds]=regexp(filename{file_nr},'[0123456789.-]*');
    for current_structure=1:nr_structures
        for parameter_nr=1:nr_parameters
            results{((file_nr - 1) * nr_structures ) + current_structure ,parameter_nr + 1 + 2}=filename{file_nr}(parameterStarts(parameter_nr):parameterEnds(parameter_nr));
        end
        results{((file_nr - 1) * nr_structures ) + current_structure , 2} = num2str(file_nr);    %save current file number
        
        %save current structure name
        filename_start_index = strfind(parameters_text{1 + current_structure}, '/');
        filename_start_index = filename_start_index(end) + 1;
        
        filename_end_index = strfind(parameters_text{1 + current_structure}, '.json');
        filename_end_index = filename_end_index(1) - 1;
        
        file_name_temporary = parameters_text{1 + current_structure};
        results{((file_nr - 1) * nr_structures ) + current_structure , 3} = file_name_temporary(filename_start_index:filename_end_index);
    end
    
end

return