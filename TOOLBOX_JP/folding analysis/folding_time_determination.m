function mean_folding_time_data = folding_time_determination(profiles, x_range, lane_times, migration_speed_min, migration_speed_max, channel, varargin)
%% calculate mean time of reaching a migration speed value
% INPUTS:
%   profiles is cell array of lane profiles
%   x_range is x_range for profiles
%   lane_times is array of folding times of each lane
%   migration_speed_min is migration distance in units of mean ladder speed relative to pocket Position from which profile data is used
%   migration_speed_max is migration distance in units of mean ladder speed relative to pocket Position up to which which profile data is used
%   channel is which gel_data image channel to use for analysis
%       Optional:
%       verbose: if 'on', plot intermediate results and mean folding times

% OUTPUT:
%   mean_folding_time = struct with fields
%       .foldingTime = mean folding times over x_range indices;
%       .migration_speed_min = min migration speed selected for folding time calculations;
%       .migration_speed_max = max migration speed selected for folding time calculations;
%       .profile_reverse_cum_integral_mono = reverse cumulative integral of each profile starting
%           from end_point to start_point, adjusted to be monotonically increasing                                                        

%% parse input
p = inputParser;
% default for number of references for background correction

addRequired(p, 'profiles');
addRequired(p, 'x_range');
addRequired(p, 'lane_times');
addRequired(p, 'migration_speed_min');
addRequired(p, 'migration_speed_max');
addRequired(p, 'channel');

% optional parameter: verbose: plot intermediate results
default_verbose = 'off';
expected_verbose = {'on', 'off'};
% check if verbose is 'on' or 'off'
addParameter(p, 'verbose', default_verbose, @(x) any(validatestring(x, expected_verbose)));
    
parse(p, profiles, x_range, lane_times, migration_speed_min, migration_speed_max, channel, varargin{:});

% print intermediate results plots
verbose_bool = strcmp(p.Results.verbose, 'on');

%% set up ranges
%modifiy lane times to include time 0 step
lane_times = [0 lane_times];

%number of lanes in gel image
num_lanes = length(profiles);

% find start end end indices in x_range corresponding to min and max migration speeds
[~, start_index] = min(abs(x_range - migration_speed_min));
[~, end_index] = min(abs(x_range - migration_speed_max));

% number of datapoints in selected range
num_profile_datapoints = end_index - start_index + 1;

%% calculate cumulative integral of each profile starting from end_point to start_point

%reverse cumulative integral of each profile starting from end_point to start_point
profile_reverse_cum_integral = zeros(num_profile_datapoints, num_lanes);

for current_lane = 1:num_lanes
    %cumulative integral of profile, starting from last value up to respective x value
    profile_reverse_cum_integral(:,current_lane) = cumtrapz( profiles{current_lane}(end_index:-1:start_index ));
    %normalize profile_reverse_cum_integral
    %profile_reverse_cum_integral(:,current_lane) = profile_reverse_cum_integral(:,current_lane) / profile_reverse_cum_integral(end,current_lane);    
    % reverse integral
    profile_reverse_cum_integral(:,current_lane) = profile_reverse_cum_integral(end:-1:1, current_lane);

end
%add time = 0 profile_reverse_cum_integral with all values = 0
profile_reverse_cum_integral = [ zeros(num_profile_datapoints, 1) profile_reverse_cum_integral] ;

% plot profiles in selected range
if verbose_bool
    figure
    for current_lane = 1:num_lanes
        plot(x_range(start_index:end_index), profiles{current_lane}(start_index:end_index));
        hold on
    end
    title('profiles');
end
pause

% plot reverse cumulative integral for each lane
if verbose_bool
    figure
    for current_lane = 1:num_lanes
        plot(x_range(start_index:end_index), profile_reverse_cum_integral(:,current_lane + 1));
        hold on
    end
    title('profile_reverse_cum_integral before monotonization over migration');
end


% calculate monotized profiles by reducing squared error between fit curve and data curve while monotonically increasing
profile_reverse_cum_integral_mono = zeros(num_profile_datapoints, num_lanes + 1);
% fit options (increase accuracy if mean folding time calculation becomes rough)
options = optimoptions('fmincon','Display','off', 'OptimalityTolerance', 10^-4);
result_vector = zeros(1,num_lanes + 1);

for current_migration_index = num_profile_datapoints - 1:-1:1
    % normalize data so that fit options can be kept constant
    current_data = profile_reverse_cum_integral(current_migration_index,:);
    max_value = max(current_data);
    current_data = current_data ./ max_value;

    % generate temporary error function for current dataset
    temporary_error_function = make_error_function(lane_times, current_data);
    lower_bound = zeros(1,num_lanes + 1);
    upper_bound = [0 ones(1,num_lanes) * Inf];
    % fit data
    result_vector = fmincon(temporary_error_function, result_vector, [], [], [], [],...
                            lower_bound, upper_bound, [], options);

    % rescale to original values to reverse previous normalization
    % calculate fit function intensity from changes in intensity
    profile_reverse_cum_integral_mono(current_migration_index, :) = cumsum(result_vector .* max_value);
end


% plot reverse cumulative integral over time for some migration speed cutoffs
% also plot monotonically increasing fit functions
trace_colors = parula(num_profile_datapoints);
if verbose_bool
    figure
    for current_migration_index = 1:20:num_profile_datapoints
        plot(lane_times, profile_reverse_cum_integral(current_migration_index, :), 'LineStyle', '-', ...
        'Color', trace_colors(current_migration_index,:) )
        hold on
        plot(lane_times, profile_reverse_cum_integral_mono(current_migration_index, :), 'LineStyle', ':', ...
        'Color', trace_colors(current_migration_index,:) )

    end
    title('profile_reverse_cum_integral');
end

    function error_function_handle = make_error_function(times, intensity_real)
        
        error_function_handle = @error_function;
        % function calculating error between a set of intensitiy increase in a monotonic function
        function error = error_function(intensity_steps_fit)

            % calculate monotonically increasing fit function
            intensities_fit = cumsum(intensity_steps_fit);

            % difference between fit and real data at lane timesteps
            diff_intensities = intensity_real - intensities_fit;
            % lengths of time intervals between lanes
            diff_times = times(2:end) - times(1:end-1);
            % number of intervals between two datapoints
            nr_intervals = length(times) - 1;

            % calculate integral over square of difference between fit curve and data curve over time interval
            error = 0;
            for step = 1:nr_intervals
                error = error + diff_times(step) * diff_intensities(step)^2 ...
                              + diff_times(step) * diff_intensities(step) * (diff_intensities(step + 1) - diff_intensities(step)) ...
                              + diff_times(step) * 1/3 * (diff_intensities(step + 1) - diff_intensities(step))^2;
            end
            %error = sqrt(error);
        end
    end

%% calculate mean folding time

%mean folding time of structures that fold until the last timestep, starting from fastest migration speed to slowest
mean_folding_time = zeros(num_profile_datapoints, 1);
for current_lane = 1 : num_lanes
    %probability of folding in current timestep
    slope = profile_reverse_cum_integral_mono(:,current_lane + 1) - profile_reverse_cum_integral_mono(:,current_lane);
    slope = slope / (lane_times(current_lane + 1) - lane_times(current_lane));
    %normalize probability distribution to 1
    slope = slope ./ profile_reverse_cum_integral_mono(:,end);
    %mean folding time contribution of timestep
    mean_folding_time = mean_folding_time + 0.5 * slope * (lane_times(current_lane + 1)^2 - lane_times(current_lane)^2);
end

if verbose_bool
    figure
    title('mean folding time');
    plot(x_range(start_index:end_index), mean_folding_time);
end

mean_folding_time_data.('foldingTime') = mean_folding_time;
mean_folding_time_data.('migration_speed_min') = migration_speed_min;
mean_folding_time_data.('migration_speed_max') = migration_speed_max;
mean_folding_time_data.('profile_reverse_cum_integral_mono') = profile_reverse_cum_integral_mono(:,end);

end