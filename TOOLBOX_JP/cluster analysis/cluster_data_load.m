function [ results ] = cluster_data_load(nr_structures,nr_time_steps, varargin)
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
nr_parameters=12;                                                   %number of parameters in filename

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
    current_line = fgetl(current_file);
    trace_nr = 0;
    
    while current_line~=-1                                         %load simulation data line by line
        for step_nr=1:1+nr_time_steps
            current_line=fgetl(current_file);
            current_data(step_nr,:)=current_data(step_nr,:)+sscanf(current_line,'%f',[1,24]);
            %step_nr
        end
        current_line = fgetl(current_file);
        trace_nr = trace_nr + 1; 
    end   

    fclose(current_file);
    
    results{file_nr,1} = current_data/trace_nr;

    %load parameter valuestart and end positions from file name
    [parameterStarts, parameterEnds]=regexp(filename{file_nr},'[0123456789.-]*');
    for parameter_nr=1:nr_parameters
        results{file_nr,parameter_nr + 1}=filename{file_nr}(parameterStarts(parameter_nr):parameterEnds(parameter_nr));
    end
    
end

return
    
%     %analyse data
%     
%     for j=0:(nr_traces-1)      %create averaged Data for each structure. empty (=0) timesteps are set to 'folded' value of 0.8
%         averagedData(:,:,file_nr)=averagedData(:,:,file_nr)+current_data(j*nr_times_steps+1:(j+1)*nr_times_steps,2:end)/nr_traces;
%         averagedData(:,:,file_nr)=averagedData(:,:,file_nr)+0.8*(current_data(j*nr_times_steps+1:(j+1)*nr_times_steps,2:end)==0)/nr_traces;
%     end
%     
%     results{file_nr+1,1}=filename{file_nr};                                     %write filenames in results cell
%                                                                     
%     [parameterStarts, parameterEnds]=regexp(results{file_nr+1,1},'[0123456789.-]*');%load parameters from file name
%     for j=1:nr_parameters
%         results{file_nr+1,j+1}=results{file_nr+1,1}(parameterStarts(j):parameterEnds(j));
%     end
%     
%     for j=1:nr_structures                              %set number of folded structures to max, avg folding time to 0, avg fold ratio to 0
%         results{file_nr+1,j+1+nr_parameters}=nr_traces;
%         results{file_nr+1,j+1+nr_parameters+nr_structures}=0;
%         results{file_nr+1,j+1+nr_parameters+2*nr_structures}=0;
%     end
%     
%     tempSelector=(current_data(:,2:end)~=0);                     %selector selects all steps with folding ratio != 0
%     thresholdSelector=(current_data(:,2:end)>=threshold);        %selector selects all steps with folding ratio above threshold
% 
%     for k=1:nr_structures                                                  
%         for j=1:nr_traces
%             if sum(current_data(j*nr_times_steps-5:j*nr_times_steps,k+1))~=0  %final timesteps not 0, simulation reached timelimit
%                 results{file_nr+1,k+1+nr_parameters}=results{file_nr+1,k+1+nr_parameters}-1;   %reduce number of folded structures if structure is unfolded at end of trace
%                 
%                 %add fold ratio at end of trace
%                 results{file_nr+1,k+1+nr_parameters+2*nr_structures}=results{file_nr+1,k+1+nr_parameters+2*nr_structures}+current_data((j-1)*nr_times_steps+nr_times_steps,k+1); 
%                 
%                 if ~threshold_bool   %add constant non-fold time to average time of folding
%                     results{file_nr+1,k+1+nr_parameters+nr_structures}=results{file_nr+1,k+1+nr_parameters+nr_structures}+current_data(nr_times_steps,1); 
%                 end
%             else	%simulation stopped before timelimit (folded/maxNrSteps)
%                 tempCurrentData=current_data((j-1)*nr_times_steps+1:j*nr_times_steps,:);     %current trace data
%                 tempCurrentData=tempCurrentData(tempSelector((j-1)*nr_times_steps+1:j*nr_times_steps,k),:);     %all timesteps where folding ratio !=0
%                 if ~isempty(tempCurrentData)        %check if structure ever started folding, add last timestep that was still folding
%                     if ~threshold_bool               %add time of folding to average time of folding
%                         results{file_nr+1,k+1+nr_parameters+nr_structures}=results{file_nr+1,k+1+nr_parameters+nr_structures}+tempCurrentData(end,1);
%                     end
%                     %add foldratio at end of folding to check simulation limits
%                     results{file_nr+1,k+1+nr_parameters+2*nr_structures}=results{file_nr+1,k+1+nr_parameters+2*nr_structures}+tempCurrentData(end,k+1);
%                 end              
%             end
%             
%             if threshold_bool
%                 tempCurrentData=current_data((j-1)*nr_times_steps+1:j*nr_times_steps,:);    %current trace data
%                 tempCurrentData=tempCurrentData(thresholdSelector((j-1)*nr_times_steps+1:j*nr_times_steps,k),:); %all timesteps where folding ratio >= threshold
%                 if ~isempty(tempCurrentData) %threshold never reached, add constant non-fold time
%                     results{file_nr+1,k+1+nr_parameters+nr_structures}=results{file_nr+1,k+1+nr_parameters+nr_structures}+tempCurrentData(1,1); 
%                 else    %add first timestep where folding ratio > threshold
%                     results{file_nr+1,k+1+nr_parameters+nr_structures}=results{file_nr+1,k+1+nr_parameters+nr_structures}+current_data(nr_times_steps,1); %add threshold crossing time
%                 end
%             end
%             
%         end
%         results{file_nr+1,k+1+nr_parameters+nr_structures}=results{file_nr+1,k+1+nr_parameters+nr_structures}/nr_traces; %calculate average time of folding
%         results{file_nr+1,k+1+nr_parameters+2*nr_structures}=results{file_nr+1,k+1+nr_parameters+2*nr_structures}/nr_traces; %calculate average folding ratio at end of simulation
%     end
%    
% end
%                                                                     %write column names from parameters text block
% results{1,1}='filename';
% 
% if ~isempty(structure_names)                                         %write structure names supplied by user
%     for i=1:nr_structures
%     results{1,i+1+nr_parameters}=structure_names{i};
%     results{1,i+1+nr_parameters+nr_structures}=structure_names{i};
%     end
% else                                                                %read names of structure from files analyzed, write to results
%     for i=1:nr_structures
%     tempLocation=strfind(parameters_text{i+1}, '/');
%     results{1,i+1+nr_parameters}=parameters_text{i+1}(tempLocation(size(tempLocation,2))+1:size(parameters_text{i+1},2));
%     results{1,i+1+nr_parameters+nr_structures}=parameters_text{i+1}(tempLocation(size(tempLocation,2))+1:size(parameters_text{i+1},2));
%     end
% end
% 
%     
% 
% %% save results to  results.txt
% 
% current_file=fopen([pathname filesep 'results.txt'],'w');
% 
% for i=1:numFiles+1
%     for j=1:2*nr_structures+1+nr_parameters
%         if ischar(results{i,j})
%             fprintf(current_file,'%s\t',results{i,j});
%         else
%             fprintf(current_file,'%f\t',results{i,j});
%         end
%     end
%     fprintf(current_file,'\n');
% end
% 
% fclose(current_file);
% 
% return
% 
% %% visualisation
% 
% parameterTypes=cell(1,numParameters);
% numParameterType=cell(1,numParameters);
% for i=1:numParameters
%     parameterTypes{i}=unique(resultParameters(2:end,i));
%     numParameterType{i}=length(parameterTypes{i});
% end
%                                                                 %calculate avg fold time of each structure for one given parameter
% temperatures=unique(results(2:end,2));
%                                                                 
% currentFigure=figure;       
% currentAxes=gca;
% for i=1:numParameters
%     tempAverages=zeros(numParameterType{i},length(temperatures));
%     for j=1:numParameterType{i}
%         for k=1:length(temperatures)
%             tempAverages(j,k)=mean([results{strcmp(resultParameters(:,i),parameterTypes{i}(j)) & strcmp(resultParameters(:,temperatureNameLocation),temperatures(k)),4}]);
%         end
%     end
%     plot(transpose(tempAverages),'LineWidth',5) %plot fold time over temperature for one constant parameter averaged over all other parameters
%     legend(parameterTypes{i})
%     ylabel('folding time [sec]')
%     xlabel('Temperature [C]')
%     axis([-inf inf 0 currentData(timeSteps,1)])
%     set(currentAxes, 'XTick',1:length(temperatures), 'XTickLabel',temperatures)
%     set(currentAxes, 'FontSize',20)
%     pause
%     
%     plot(tempAverages,'LineWidth',5) %plot fold time over one parameter for constant temperature averaged over all other parameters
%     legend(temperatures)
%     ylabel('folding time [sec]')
%     xlabel('parameter')
%     axis([-inf inf 0 currentData(timeSteps,1)])
%     set(currentAxes, 'XTick',1:numParameterType{i}, 'XTickLabel',parameterTypes{i})
%     set(currentAxes, 'FontSize',20)
%     pause
% end
% 
% num=1:200;                                                       %parameter file counter
% 
% currentFigure=figure;                                               %plot all results
% currentBar=bar3(cell2mat(results(2:end,2+numStructures+1:2+2*numStructures)));
% for k = 1:length(currentBar)                                        %set color value to bar height
%     zdata = get(currentBar(k),'ZData');
%     set(currentBar(k),'CData',zdata,...
%         'FaceColor','interp')
% end
% set(gca,'CLimMode','auto')                                          %set color range to auto full range
% view(0, 90);                                                        %rotate to top view
% title('All');
% set(gca, 'XTickLabel',structureNames)
% xlabel('structure')
% ylabel('parameter Set')
% colorbar
% figAll=currentFigure;
% 
% 
% currentFigure=figure;                                               %plot all temperatures
% currentBar=bar3(cell2mat(results(2:end,2+numStructures+1:2+2*numStructures)));
% for k = 1:length(currentBar)                                        %set color value to bar height
%     zdata = get(currentBar(k),'ZData');
%     set(currentBar(k),'CData',zdata,...
%         'FaceColor','interp')
% end
% currentAxes=get(currentFigure,'CurrentAxes');
% set(currentAxes,'CLimMode','auto')                                  %set color range to auto full range
% view(0, 90);                                                        %rotate to top view
% caxis([0 7200])
% title('one set');                                                   %set title
% set(currentAxes, 'XTickLabel',structureNames)                       %set x axis labels to structure names
% set(currentAxes, 'YTickLabel',results(2:end,2))                     %set y axis labels to temperature values
% xlabel('structure')                                                 %label x axis
% ylabel('temperature [C]')                                          %label y axis
% cb = colorbar('vert');                                              %add a colormap legend
% zlab = get(cb,'ylabel');                                            %set colormap legend label
% set(zlab,'String','Average Fold Time [sec]'); 
% colormap(flipud(copper))                                            %set colormap
% set(findall(currentFigure,'-property','FontSize'),'FontSize',18)    %change all font sizes
% 
% set(currentFigure, 'Units', 'centimeters')                          %set paper size to figure size
% set(currentFigure, 'PaperUnits', 'centimeters')
% figureSize=get(currentFigure, 'Position');
% set(currentFigure, 'PaperSize', figureSize(3:4))
% 
% 
% set(0,'defaultAxesLineStyleOrder','-|--|:')
% currentFigure=figure;                                               %line plot folding times
% plot(cell2mat(results(2:end,2+numStructures+1:2+2*numStructures)),'LineWidth',3)
% currentAxes=get(currentFigure,'CurrentAxes');
% axis([-inf inf 0 7200])
% set(currentAxes, 'XTick',1:length(results(2:end,2)), 'XTickLabel',cell2mat(results(2:end,2)))
% legend(structureNames);
% xlabel('temperature [C]')                                           %label x axis
% ylabel('Average Fold Time [sec]')                                          %label y axis
% set(findall(currentFigure,'-property','FontSize'),'FontSize',18)    %change all font sizes
% 
% set(currentFigure, 'Units', 'centimeters')                                    %set paper size to figure size
% set(currentFigure, 'PaperUnits', 'centimeters')
% figureSize=get(currentFigure, 'Position');
% set(currentFigure, 'PaperSize', figureSize(3:4))
% 
% 
% %resultsTempSort=zeros(0,0)
% 
% temperatures=unique(results(2:end,2));
% for i=1:length(temperatures);                                       %plot results for each temperature
%     
%     %resultsTempSort=[resultsTempSort cell2mat(results(strcmp(results(:,2),temperatures(i)),2+numStructures+1:2+2*numStructures))];
%     
%     currentFigure=figure;
%     currentBar=bar3(cell2mat(results(strcmp(results(:,2),temperatures(i)),2+numStructures+1:2+2*numStructures)));
%     for k = 1:length(currentBar)
%         zdata = get(currentBar(k),'ZData');
%         set(currentBar(k),'CData',zdata,...
%             'FaceColor','interp')
%     end
%     set(gca,'CLimMode','auto')
%     view(0, 90);
%     title(temperatures(i));
%     set(gca, 'XTickLabel',structureNames)
%     set(gca, 'YTickLabel',num(strcmp(results(:,2),temperatures(i))))
%     xlabel('structure')
%     ylabel('parameter Set')
%     colorbar
% end