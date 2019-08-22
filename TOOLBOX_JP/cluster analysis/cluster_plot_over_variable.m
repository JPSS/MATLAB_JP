function [ output_args ] = cluster_plot_over_variable( data, variable_index, plot_y_max, variable_in_one_plot, variable_in_one_curve, varargin)
%% plot average traces of structures with same parameters at different temperatures
%   data is results from cluster_data_load()
%   variable_index is index of simulation variable to plot
%   plot_y_max is y plotrange
%   variable_in_one_plot is variable that is held constant in a single subplot
%   variable_in_one_curve is variable that is held constant in a single curve

%variable_index=19;  %bases good
%variable_index=20;  %bases bad

%variable_index=8; %crossover unique ratio
%variable_index=16; %crossovers ratio
%variable_index=13; %crossovers good
%variable_index=14; %crossovers good ratio
%variable_index=21; %crossovers good unique ratio

%variable_index=17; %number bound staple
%variable_index=22; %number basestacks good
%variable_index=23; %number segments good
%variable_index=24; %phi critical

%% parse input variables
parser = inputParser;

% required parameter
addRequired(parser,'data');
addRequired(parser,'variable_index');
addRequired(parser,'plot_y_max');
addRequired(parser,'variable_in_one_plot');
addRequired(parser,'variable_in_one_curve');

% optional parameter: 
addParameter(parser,'timesteps',{});

parse(parser, data, variable_index, plot_y_max, variable_in_one_plot, variable_in_one_curve, varargin{:});

timesteps = parser.Results.timesteps

%% determine changing parameters

number_of_parameters = size(data,2);
parameters_in_parameter_set = [3:number_of_parameters];
parameters_in_parameter_set = parameters_in_parameter_set(find(parameters_in_parameter_set ~= variable_in_one_curve));
parameters_in_parameter_set = parameters_in_parameter_set(find(parameters_in_parameter_set ~= variable_in_one_plot));

%determine different values for each parameter
parameter_values = {};
number_of_parameter_values = [];
for parameter = 2 : number_of_parameters
    parameter_values{parameter} = unique(data(:,parameter));
    [~, sorting_indices] = sort(str2double(parameter_values{parameter}));
    parameter_values{parameter} = parameter_values{parameter}(sorting_indices);
    number_of_parameter_values(parameter) = length(parameter_values{parameter});
end
    
%parameter indices of current parameter set, inculidng unchanging parameters
current_parameter_indices = ones(1,number_of_parameters);
%current changing parameter, counting includes unchanging parameters
changing_parameter = 1;


%% plot data

figure('pos',[0,0,600,400])

%cycle through all parameter sets
while changing_parameter <= length(parameters_in_parameter_set)
    
    %select all results with same parameter set
    subset_selector = ones(1,size(data,1));
    for parameter = parameters_in_parameter_set
        subset_selector = subset_selector & ismember( {data{:,parameter}} , { parameter_values{parameter}{current_parameter_indices(parameter)} });
    end

    %print current parameter set
    results_subset = data(subset_selector,:);
    if ~isempty(results_subset)
        strcat({'param Nr = ','exp = ','mult = ','overhang = ','stack = ','dangle = ','temp = ','minLen = '},results_subset(1 , [2,5,7:12]) )
    end
    
    %current subplot number in subplot plot
    subplot_index=1;
    clf
    
    % total number of subplots
    number_of_subplots = length(parameter_values{variable_in_one_plot});
    % number of subplots in x or y direction
    subplots_per_axis = ceil(sqrt(number_of_subplots));
    
    %select subset to plot from parameter set subset
    for plot_parameter_value = parameter_values{variable_in_one_plot}'
        subset_to_plot = subset_selector & ismember( {data{:,variable_in_one_plot}} , plot_parameter_value );
        results_subset = sortrows( data(subset_to_plot,:), variable_in_one_curve);

        subplot(subplots_per_axis, subplots_per_axis, subplot_index)
        plot_x_range_current = 0;

        %plot data from current subset
        for current_result_index = 1 : size(results_subset,1)
            set(gca,'LineStyleOrder','-|--|:')
            plot( results_subset{ current_result_index, 1 }( :, variable_index ), 'DisplayName' , num2str(results_subset{ current_result_index, variable_in_one_curve }) )
            hold on
            plot_x_range_current = adjust_range( results_subset{ current_result_index }, plot_x_range_current );
            %set title to overhang penalty, stacking, multiplier, temperature, min len
            title( strcat({'mult = ','overhang = ','stack = ','dangle = ','temp = ','minLen = '}, results_subset(1,[7:12]) ) )
        end
        legend({results_subset{:, 3}});
        subplot_index = subplot_index + 1;
    end
    pause
    
    %cycle to next parameter set
    increased_parameter = 0;
    while ~increased_parameter
        
        %check if current parameter index can be increased
        if current_parameter_indices( parameters_in_parameter_set(changing_parameter) ) < number_of_parameter_values( parameters_in_parameter_set(changing_parameter) )
            %increase current parameter index
            current_parameter_indices( parameters_in_parameter_set(changing_parameter) ) = current_parameter_indices( parameters_in_parameter_set(changing_parameter) ) + 1;
            %set previous parameter indices to 1 again
            for parameter = 1 : changing_parameter - 1
                current_parameter_indices( parameters_in_parameter_set(parameter) ) = 1;
            end
            %set changing parameter to first again
            changing_parameter = 1;
            %exit loop
            increased_parameter = 1;
        else
            %switch to next parameter if previous parameter index can't be increased
            changing_parameter = changing_parameter + 1;
        end
        
        %abort if changing parameter goes beyond number of parameters in parameter set
        if changing_parameter > length(parameters_in_parameter_set)
            increased_parameter = 1;
        end
    end
end
        
%determine plot range for current trace
function range_new = adjust_range(current_result,range_current)
    %find first 0 in plotted variable after 50th value
    range_new = find(current_result( 50 : end, variable_index ) == 0, 1, 'first');
    if ~isempty(range_new)
        range_new = range_new + 50;
        %if timesteps parameter was passed, set fixed x axis lenght
        if ~isempty(timesteps)
            range_new = timesteps;
        end
        if range_new > range_current
            axis([0 range_new 0 plot_y_max])
        else
            range_new = range_current;
        end
    else
        range_new = length(current_result(:, variable_index ));
        axis([0 range_new 0 plot_y_max])
    end
end

end