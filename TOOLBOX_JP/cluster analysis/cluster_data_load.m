function [ results , mean_times ] = cluster_data_load(nr_structures,nr_time_steps, varargin)
%% analyzes gillespie simulation data
%   nr_structures is number of structures in file
%   nr_times_steps is number of steps printed

%   filename{N} are filenames of files
%   pathname is pathnames of files

%% parse input variables
    p = inputParser;
    % required parameter
    addRequired(p,'nr_structures');
    addRequired(p,'nr_times_steps');
    
    
%nr of steps to be considered in mean folding time calculation
nr_mean_folding_time_steps = 10;
    
%% load data
%opens file, copies data into results and currentData, analyzes,
%then opens next file
%current_data is array of averaged traces from file

[filename, pathname]=uigetfile('*','select cluster data','MultiSelect','on');

nr_parameter_lines = 28;                                               %number of lines after structure names before data starts
parameters_text = cell(1+nr_structures+nr_parameter_lines,1);           %saves parameters text, each cell one line
nr_parameters = 12;                                                   %number of parameters in filename

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
   
    current_data=zeros(nr_time_steps+1,24);
    %current mean folding times
    current_mean_times = zeros(1, nr_mean_folding_time_steps);
    
    current_line = fgetl(current_file);
    trace_nr = 0;
    
    %load simulation data line by line
    while current_line~=-1
        %current folding times
        current_times = nr_time_steps * ones(1, nr_mean_folding_time_steps);
        for step_nr=1:1+nr_time_steps
            current_line=fgetl(current_file);
            current_line_data = sscanf(current_line,'%f',[1,24]);
            current_data(step_nr,:)=current_data(step_nr,:) + current_line_data;
            
            for cutoff = 1:nr_mean_folding_time_steps
                if current_line_data(14) > cutoff/nr_mean_folding_time_steps && current_times(cutoff) == nr_time_steps
                    current_times(cutoff) = step_nr;
                end
            end
        end
        current_mean_times = current_mean_times + current_times;
        current_line = fgetl(current_file);
        trace_nr = trace_nr + 1; 
    end   

    fclose(current_file);
    
    results{file_nr,1} = current_data/trace_nr;
    mean_times(file_nr,:) = current_mean_times / trace_nr;

    %load parameter valuestart and end positions from file name
    [parameterStarts, parameterEnds] = regexp(filename{file_nr},'[0123456789.-]*');
    for parameter_nr = 1:nr_parameters
        results{file_nr, parameter_nr + 1} = filename{file_nr}(parameterStarts(parameter_nr):parameterEnds(parameter_nr));
    end
    
    %save current structure name
    current_structure_number = str2num( results{file_nr, 3} ) + 1;
    filename_start_index = strfind(parameters_text{1 + current_structure_number}, '/');
    filename_start_index = filename_start_index(end) + 1;

    filename_end_index = strfind(parameters_text{1 + current_structure_number}, '.json');
    filename_end_index = filename_end_index(1) - 1;

    file_name_temporary = parameters_text{1 + current_structure_number};
    results{file_nr, 3} = file_name_temporary(filename_start_index : filename_end_index);    
    
end

return