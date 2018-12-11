function gel_data = folding_time_determination(gel_data, lane_times, migration_speed_min, migration_speed_max, num_points, channel, varargin)
%% Loads image, fits lanes according to step function convolved with gaussian
%   INPUTS:
%       gel_data struct after applying analyze_Mean_ladder() to it
%       lane_times is array of folding times of each lane
%       migration_speed_min is migration distance in units of mean ladder speed relative to pocket Position from which profile data is used
%       migration_speed_max is migration distance in units of mean ladder speed relative to pocket Position up to which which profile data is used
%       num_points is number of datapoints in rescaled profiles
%       channel is which gel_data image channel to use for analysis

%   OUTPUT:
%
%% parse input
p = inputParser;
% default for number of references for background correction

addRequired(p, 'gel_data');
addRequired(p, 'lane_times');
addRequired(p, 'migration_speed_min');
addRequired(p, 'migration_speed_max');
addRequired(p, 'num_points');
addRequired(p, 'channel');

% optional parameter: automatic fit starting value selection
default_verbose = 'off';
expected_verbose = {'on', 'off'};
% check automatic_start_values is 'on' or 'off'
addParameter(p, 'verbose', default_verbose, @(x) any(validatestring(x,expected_verbose)));
    
parse(p, gel_data, lane_times, migration_speed_min, migration_speed_max, num_points, channel, varargin{:});

% print intermediate results plots
verbose_bool = strcmp(p.Results.verbose, 'on');

%% set up ranges
%modifiy lane times to include time 0 step
lane_times = [0 lane_times];

%number of lanes in gel image
num_lanes = length(gel_data.profiles);
% number of datapoints in original profiles
num_profile_datapoints = length(gel_data.profiles{1,1});

mean_ladder_speeds = gel_data.ladder_correction.mean_ladder_speeds;

%new data range from and including start_point to and including end_point in units of mean ladder speed
new_data_range = linspace(migration_speed_min, migration_speed_max, num_points);
%profiles rescaled onto new data range
rescaled_profiles = zeros(num_points, num_lanes);


%% calculate cumulative integral of each profile starting from end_point to start_point

%reverse cumulative integral of each profile starting from end_point to start_point
profile_reverse_cum_integral = zeros(num_points, num_lanes);

for i = 1:num_lanes
    
    %x values of profile data in units of corrected mean ladder speed
    old_data_range = ((1:num_profile_datapoints) - gel_data.pocketPositions(i)) / mean_ladder_speeds(i);
    
    %interpolated profile data in selected subrange in mean ladder speed units
    rescaled_profiles(:,i) = interp1(old_data_range, gel_data.profiles{channel, i}, new_data_range);
    
    %multiply new profiles with ladder speed to conserve total intensity integral value 
    rescaled_profiles(:,i) = rescaled_profiles(:,i) * mean_ladder_speeds(i);
    
    %cumulative integral of profile, starting from last value
    profile_reverse_cum_integral(:,i) = cumtrapz(fliplr(rescaled_profiles(end:-1:1,i)));
    %normalize profile_reverse_cum_integral
    profile_reverse_cum_integral(:,i) = profile_reverse_cum_integral(:,i) / profile_reverse_cum_integral(end,i);
end

if verbose_bool
    figure
    title('rescaled_profiles');
    for current_lane = 1:num_lanes
        plot(linspace(migration_speed_max,migration_speed_min,num_points),rescaled_profiles(end:-1:1,current_lane));
        hold on
    end
end


if verbose_bool
    figure
    title('profile_reverse_cum_integral before monotonization');
    for current_lane = 1:num_lanes
        plot(linspace(migration_speed_max,migration_speed_min,num_points),profile_reverse_cum_integral(:,current_lane));
        hold on
    end
end

%make profile_reverse_cum_integral monotonic increasing over time
for current_lane = 2:num_lanes
    %replace all values that are smaller than in previous lane profile_reverse_cum_integral
    indices_of_smaller_values = profile_reverse_cum_integral(:,current_lane - 1) > profile_reverse_cum_integral(:,current_lane);
    profile_reverse_cum_integral(indices_of_smaller_values,current_lane) = profile_reverse_cum_integral(indices_of_smaller_values, current_lane - 1);
end

%add time = 0 profile_reverse_cum_integral with all values=0 except for last value=1
profile_reverse_cum_integral = [ zeros(num_points,1) profile_reverse_cum_integral] ;
profile_reverse_cum_integral(end,1) = 1;

if verbose_bool
    figure
    title('profile_reverse_cum_integral after monotonization');
    for current_lane = 2:num_lanes+1
        plot(linspace(migration_speed_max,migration_speed_min,num_points),profile_reverse_cum_integral(:,current_lane));
        hold on
    end
end

%% calculate mean folding time

%mean folding time of structures that fold until the last timestep, starting from fastest migration speed to slowest
mean_folding_time = zeros(num_points,1);
for current_lane = 1 : num_lanes
    %probability of folding in current timestep
    slope = profile_reverse_cum_integral(:,current_lane + 1) - profile_reverse_cum_integral(:,current_lane);
    slope = slope / (lane_times(current_lane + 1) - lane_times(current_lane));
    %normalize probability distribution to 1
    slope = slope ./ profile_reverse_cum_integral(:,end);
    %mean folding time contribution of timestep
    mean_folding_time = mean_folding_time + 0.5 * slope * (lane_times(current_lane + 1)^2 - lane_times(current_lane)^2);
end

if verbose_bool
    figure
    title('mean folding time');
    plot(linspace(migration_speed_max,migration_speed_min,num_points),mean_folding_time);
end

gel_data.('foldingTime') = mean_folding_time;
gel_data.('migration_speed_min') = migration_speed_min;
gel_data.('migration_speed_max') = migration_speed_max;
gel_data.('profile_reverse_cum_integral_final') = profile_reverse_cum_integral(:,end);

end