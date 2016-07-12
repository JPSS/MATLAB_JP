%% this generates a cell array with the parameter sets of a given simulation run, and their respective parameter set index

array{1} = [5];                        %loop rate exponents - arrayDuplicats
array{2} = [-2.44];                    %constant dG for loop rate factor - arrayExponents
array{3} = [-0.002579978];             %constant dG for non unique (double) good crossover formation - arrayDGLoop
array{4} = [0.0077, 0.0088, 0.01, 0.0111, 0.0122];%multiplicator double crossover premium - arrayMultiplicator
array{5} = [12.17222222, 16.15921667]; %energy penalty touching staple overhangs - arrayDGOverlap
array{6} = [-8.33, -14.52222222];    %stacking energy - arrayDGBasestack
array{7} = [-2.534, 0];                 %dangling ends energy - arrayDGDangle
array{8} = [50, 55];                %temperature - arrayTemperatures
array{9} = [20, 7, 4];                 %minimum random match length - arrayMinLength
array{10} = [10000];                   %update step for forced loop rate update - arrayUpdateStep

nr_parameters = length(array);

%current position in each parameter array
current_positions = ones(nr_parameters,1);
current_parameter_combination = 0;

%output result
results={};

running = 1;
while running
    %print parameter values to results
    results{current_parameter_combination + 1,1} = current_parameter_combination;
    for level = 1:nr_parameters
        results{current_parameter_combination + 1, level + 1} = array{ level }( current_positions( level ) );
    end
    
    %advance index parameter combination
    current_parameter_combination = current_parameter_combination + 1;
    %advance topmost parameter position by 1
    current_positions(nr_parameters) = current_positions(nr_parameters) + 1;
    

    %advance current parameter positions
    level = nr_parameters;
    while current_positions(level) == length(array{level}) + 1
        
        %reached end of current parameter array, set to position to 1 again
        current_positions(level) = 1;
        %check level below
        level = level - 1;
        
         %check if reached final parameter set
        if level == 0;
            running = 0;
            break
        end
        
        %increase parameter set below by 1
        current_positions(level) = current_positions(level) + 1;

    end
    
end
        
    
    