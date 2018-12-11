function [ output_args ] = cluster_plot_folding_times( data, mean_times, varargin)
%% plot average traces of structures with same parameters at different temperatures
%   data is results from cluster_data_load()
%   mean_times is mean_time from cluster_data_load()
%   variable_index is index of simulation variable to plot
%   plot_y_max is y plotrange
%   time_index is column from mean_times to plot
%   variable_in_one_plot is constant variable in one plot
%   nr_plots_x is number of subplots in x direction
%   nr_plots_y is number of subplots in y direction
%   plotstyle (single: new plot for each parameter set, joint: all plots in one graph)
%   names is cell of strings, names of structures in simulation
%   logarithm is string on or off logarithm time scale

%% parse input variables
parser = inputParser;
% required parameter
addRequired(parser,'data');
addRequired(parser,'mean_times');

% optional parameter: 
addParameter(parser,'plot_y_min',0);
addParameter(parser,'plot_y_max',0);
addParameter(parser,'time_index',5);
addParameter(parser,'variable_in_one_plot',11);
addParameter(parser,'nr_plots_x',4);
addParameter(parser,'nr_plots_y',3);
addParameter(parser,'names',{});

default_plotstyle = 'single';
expected_plotstyle = {'single', 'joint'};
addParameter(parser,'plotstyle', default_plotstyle,  @(x) any(validatestring(x,expected_plotstyle)));

default_logarithmic = 'off';
expected_logarithmic = {'on', 'off'};
addParameter(parser,'logarithmic', default_logarithmic,  @(x) any(validatestring(x,expected_logarithmic)));

default_legend = 'on';
expected_legend = {'on', 'off'};
addParameter(parser,'legend', default_legend,  @(x) any(validatestring(x,expected_legend)));

default_skip_empty = 'off';
expected_skip_empty = {'on', 'off'};
addParameter(parser,'skip_empty', default_skip_empty,  @(x) any(validatestring(x,expected_skip_empty)));

parse(parser, data, mean_times, varargin{:});

plot_y_min = parser.Results.plot_y_min;
plot_y_max = parser.Results.plot_y_max;
time_index = parser.Results.time_index;
variable_in_one_plot = parser.Results.variable_in_one_plot;
nr_plots_x = parser.Results.nr_plots_x;
nr_plots_y = parser.Results.nr_plots_y;
names = parser.Results.names;

parser.Results.plotstyle;
parser.Results.logarithmic;

%% determine changing parameters

number_of_parameters = size(data,2)
parameters_in_parameter_set = [4:number_of_parameters];
parameters_in_parameter_set = parameters_in_parameter_set(find(parameters_in_parameter_set ~= variable_in_one_plot));
%variable that is held constant in a single curve
variable_in_one_curve = 3;

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
figure('pos',[0,0,900,900])
%determine numerical temperature values for current parameter set
set_variable_range = sort( cellfun( @str2num , parameter_values{variable_in_one_plot} ) );

parameter_set_nr = 1;

%cycle through all parameter sets
while changing_parameter <= length(parameters_in_parameter_set)
    
    if strcmp(parser.Results.plotstyle, 'single')
        %clear figure
        clf
    elseif strcmp(parser.Results.plotstyle, 'joint')
        subplot(nr_plots_y, nr_plots_x, parameter_set_nr)
    end
  
    %select all results with same parameter set
    subset_selector = ones(1,size(data,1));
    for parameter = parameters_in_parameter_set
        subset_selector = subset_selector & ismember( {data{:,parameter}} , { parameter_values{parameter}{current_parameter_indices(parameter)} });
    end
    
    %select plot range
    %fastest mean time
    smallest_mean_time = min(mean_times(subset_selector,time_index));
    %slowest mean time
    largest_mean_time = max(mean_times(subset_selector,time_index));
    %check if there are any traces in current parameter set
    if isempty(smallest_mean_time)
        smallest_mean_time = 1;
    end
    if isempty(largest_mean_time)
        largest_mean_time = Inf;
    end

    %print current parameter set
    results_subset = data(subset_selector,:);
    if ~isempty(results_subset)
        strcat({'param Nr = ','exp = ','mult = ','overhang = ','stack = ','dangle = ','minLen = '},results_subset(1 , [2,5,7:10,12]) )
    end

    %select subset to plot from parameter set subset
    for plot_parameter_value = parameter_values{variable_in_one_curve}'
        %select subset to plot in once curve
        subset_to_plot = subset_selector & ismember( {data{:,variable_in_one_curve}} , plot_parameter_value );
        %sort subset by variable in plot
        [results_subset, sorting_indices ] = sortrows( data(subset_to_plot,:), variable_in_one_plot);
        %select subset of mean_times and sort them as well
        sorted_times = mean_times(subset_to_plot,:);
        sorted_times = sorted_times(sorting_indices,:);

        set(gca,'LineStyleOrder','-o|--o|:o')
        
        %set y scale to logarithmic
        if strcmp(parser.Results.logarithmic, 'on')
            set(gca, 'YScale', 'log')
        end
        
        hold on
        plot( cellfun( @str2num,results_subset(:,variable_in_one_plot) ), sorted_times(:,time_index), 'DisplayName' , num2str(plot_parameter_value{1}) )
            colors = jet(10);
      
        %set y axis range
        edge_with = (max(set_variable_range) - min(set_variable_range)) / 10;
        % if variable range is size 0, set a minimum edge with
        if edge_with == 0
            edge_with = 0.1 * max(set_variable_range);
        end
        
        if ~plot_y_max
            axis([min(set_variable_range) - edge_with,  max(set_variable_range) + edge_with,...
                max(0.01, plot_y_min), max(plot_y_min+0.1, 10 * smallest_mean_time)])
        else
            axis([min(set_variable_range) - edge_with,  max(set_variable_range) + edge_with,...
                max(0.01, plot_y_min), plot_y_max])
        end
        if strcmp(parser.Results.logarithmic, 'on')
            axis([min(set_variable_range) - edge_with,  max(set_variable_range) + edge_with,...
                smallest_mean_time * 0.99, largest_mean_time * 1.01])
        end
        
        hold on
        if ~isempty(results_subset)
            title( strcat({'mult = ','overhang = ','stack = ','dangle = ','minLen = '}, results_subset(1,[7:10,12]) ) )
        end

    end
    if strcmp(parser.Results.legend, 'on')
        legend(names, 'Interpreter', 'none')
    end
    
    if strcmp(parser.Results.plotstyle, 'single')
        pause
    end
    
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
    
    if strcmp(parser.Results.skip_empty, 'off') || ~isempty(results_subset)
        parameter_set_nr = parameter_set_nr + 1;
    end
end

end