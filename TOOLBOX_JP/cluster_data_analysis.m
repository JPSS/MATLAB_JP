function [ results, averagedData ] = cluster_data_analysis(numStructures,traces,timeSteps, varargin)
%% analyzes gillespie simulation data
%   numStructures is number of structures simulated
%   traces is number of traces per structure and parameter set
%   timeSteps is number of steps printed
%   varargin is cell of supplied structure names

%   filename{N} are filenames of files
%   pathname is pathnames of files

[filename, pathname]=uigetfile('*','select cluster data','MultiSelect','on');

if ~isempty(varargin)                                                            % names of structures
    structureNames=varargin{1};
end

%% load data and analyze
%opens file, copies data into results and currentData, analyzes,
%then opens next file
%currentData is array of traceData x structures, with empty lines between
%traces omited

numParameterLines=27;                                               %number of lines after structure names before data starts
parametersText=cell(1+numStructures+numParameterLines,1);           %saves parameters text, each cell one line
numParameters=10;                                                   %number of parameters in filename

numFiles=size(filename,2);                                          %number of files loaded

if ~(iscell(filename))                                              %check if a single file was loaded, if yes, adjust data type of filename
    numFiles=1;
    filename={filename};
end

currentData=zeros(timeSteps*traces,numStructures+1);                %temporary storage current file
results=cell(numFiles+1,1+numParameters+2*numStructures);                         %analysis results of data
averagedData=zeros(timeSteps,numStructures,numFiles);

for i=1:numFiles
    
    currentFile=fopen([pathname filename{i}],'r');                  %load current file
    filename{i}
    for j=1:1+numStructures+numParameterLines
        parametersText{j}=fgetl(currentFile);                           %load parameter text
    end

    for j=0:(traces-1)                                              %load simulation data
        fgetl(currentFile);                             
        for k=1:timeSteps
            currentData(k+j*timeSteps,:)=fscanf(currentFile,'%f',numStructures+1);
        end
    end   

    fclose(currentFile);
    
    %analyse data
    for j=0:(traces-1)
        averagedData(:,:,i)=averagedData(:,:,i)+currentData(j*timeSteps+1:(j+1)*timeSteps,2:end)/traces;
        averagedData(:,:,i)=averagedData(:,:,i)+0.8*(currentData(j*timeSteps+1:(j+1)*timeSteps,2:end)==0)/traces;
    end
    
    results{i+1,1}=filename{i};                                     %put filenames in results
                                                                    
    [parameterStarts, parameterEnds]=regexp(results{i+1,1},'[0123456789.-]*');%load parameters from file name
    for j=1:numParameters
        results{i+1,j+1}=results{i+1,1}(parameterStarts(j):parameterEnds(j));
    end
    
    for j=1:numStructures                                           %set number of folded structures to max, avg folding time to 0
        results{i+1,j+1+numParameters}=traces;
        results{i+1,j+1+numParameters+numStructures}=0;
    end
    
    for k=1:numStructures                                                  
        for j=1:traces

            if sum(currentData(j*timeSteps-5:j*timeSteps,k+1))~=0  %reduce number of folded structures if structure is unfolded at end of trace
                results{i+1,k+1+numParameters}=results{i+1,k+1+numParameters}-1;                %add constant non-fold time to average time of folding
                results{i+1,k+1+numParameters+numStructures}=results{i+1,k+1+numParameters+numStructures}+currentData(timeSteps,1);
            else                                                    
                tempCurrentData=currentData((currentData((j-1)*timeSteps+1:j*timeSteps,1+k)~=0),1);     %find time of folding
                if ~isempty(tempCurrentData)                                                            %check if structure ever started folding
                    results{i+1,k+1+numParameters+numStructures}=results{i+1,k+1+numParameters+numStructures}+tempCurrentData(end); %add time of folding to average time of folding
                end              
            end
        end
        results{i+1,k+1+numParameters+numStructures}=results{i+1,k+1+numParameters+numStructures}/traces; %calculate average time of folding
    end
   
end
                                                                    %write column names from parameters text block
results{1,1}='filename';

if exist('structureNames','var')
    for i=1:numStructures                                               %read names of structure files analyzed, write to results
    results{1,i+1+numParameters}=structureNames{i};
    results{1,i+1+numParameters+numStructures}=structureNames{i};
    end
else
    for i=1:numStructures                                               %read names of structure files analyzed, write to results
    tempLocation=strfind(parametersText{i+1}, '/');
    results{1,i+1+numParameters}=parametersText{i+1}(tempLocation(size(tempLocation,2))+1:size(parametersText{i+1},2));
    results{1,i+1+numParameters+numStructures}=parametersText{i+1}(tempLocation(size(tempLocation,2))+1:size(parametersText{i+1},2));
    end
end

    

%% save results to  results.txt

currentFile=fopen([pathname filesep 'results.txt'],'w');

for i=1:numFiles+1
    for j=1:2*numStructures+1+numParameters
        if ischar(results{i,j})
            fprintf(currentFile,'%s\t',results{i,j});
        else
            fprintf(currentFile,'%f\t',results{i,j});
        end
    end
    fprintf(currentFile,'\n');
end

fclose(currentFile);

return

%% visualisation

parameterTypes=cell(1,numParameters);
numParameterType=cell(1,numParameters);
for i=1:numParameters
    parameterTypes{i}=unique(resultParameters(2:end,i));
    numParameterType{i}=length(parameterTypes{i});
end
                                                                %calculate avg fold time of each structure for one given parameter
temperatures=unique(results(2:end,2));
                                                                
currentFigure=figure;       
currentAxes=gca;
for i=1:numParameters
    tempAverages=zeros(numParameterType{i},length(temperatures));
    for j=1:numParameterType{i}
        for k=1:length(temperatures)
            tempAverages(j,k)=mean([results{strcmp(resultParameters(:,i),parameterTypes{i}(j)) & strcmp(resultParameters(:,temperatureNameLocation),temperatures(k)),4}]);
        end
    end
    plot(transpose(tempAverages),'LineWidth',5) %plot fold time over temperature for one constant parameter averaged over all other parameters
    legend(parameterTypes{i})
    ylabel('folding time [sec]')
    xlabel('Temperature [C]')
    axis([-inf inf 0 currentData(timeSteps,1)])
    set(currentAxes, 'XTick',1:length(temperatures), 'XTickLabel',temperatures)
    set(currentAxes, 'FontSize',20)
    pause
    
    plot(tempAverages,'LineWidth',5) %plot fold time over one parameter for constant temperature averaged over all other parameters
    legend(temperatures)
    ylabel('folding time [sec]')
    xlabel('parameter')
    axis([-inf inf 0 currentData(timeSteps,1)])
    set(currentAxes, 'XTick',1:numParameterType{i}, 'XTickLabel',parameterTypes{i})
    set(currentAxes, 'FontSize',20)
    pause
end

num=1:200;                                                       %parameter file counter

currentFigure=figure;                                               %plot all results
currentBar=bar3(cell2mat(results(2:end,2+numStructures+1:2+2*numStructures)));
for k = 1:length(currentBar)                                        %set color value to bar height
    zdata = get(currentBar(k),'ZData');
    set(currentBar(k),'CData',zdata,...
        'FaceColor','interp')
end
set(gca,'CLimMode','auto')                                          %set color range to auto full range
view(0, 90);                                                        %rotate to top view
title('All');
set(gca, 'XTickLabel',structureNames)
xlabel('structure')
ylabel('parameter Set')
colorbar
figAll=currentFigure;


currentFigure=figure;                                               %plot all temperatures
currentBar=bar3(cell2mat(results(2:end,2+numStructures+1:2+2*numStructures)));
for k = 1:length(currentBar)                                        %set color value to bar height
    zdata = get(currentBar(k),'ZData');
    set(currentBar(k),'CData',zdata,...
        'FaceColor','interp')
end
currentAxes=get(currentFigure,'CurrentAxes');
set(currentAxes,'CLimMode','auto')                                  %set color range to auto full range
view(0, 90);                                                        %rotate to top view
caxis([0 7200])
title('one set');                                                   %set title
set(currentAxes, 'XTickLabel',structureNames)                       %set x axis labels to structure names
set(currentAxes, 'YTickLabel',results(2:end,2))                     %set y axis labels to temperature values
xlabel('structure')                                                 %label x axis
ylabel('temperature [C]')                                          %label y axis
cb = colorbar('vert');                                              %add a colormap legend
zlab = get(cb,'ylabel');                                            %set colormap legend label
set(zlab,'String','Average Fold Time [sec]'); 
colormap(flipud(copper))                                            %set colormap
set(findall(currentFigure,'-property','FontSize'),'FontSize',18)    %change all font sizes

set(currentFigure, 'Units', 'centimeters')                          %set paper size to figure size
set(currentFigure, 'PaperUnits', 'centimeters')
figureSize=get(currentFigure, 'Position');
set(currentFigure, 'PaperSize', figureSize(3:4))


set(0,'defaultAxesLineStyleOrder','-|--|:')
currentFigure=figure;                                               %line plot folding times
plot(cell2mat(results(2:end,2+numStructures+1:2+2*numStructures)),'LineWidth',3)
currentAxes=get(currentFigure,'CurrentAxes');
axis([-inf inf 0 7200])
set(currentAxes, 'XTick',1:length(results(2:end,2)), 'XTickLabel',cell2mat(results(2:end,2)))
legend(structureNames);
xlabel('temperature [C]')                                           %label x axis
ylabel('Average Fold Time [sec]')                                          %label y axis
set(findall(currentFigure,'-property','FontSize'),'FontSize',18)    %change all font sizes

set(currentFigure, 'Units', 'centimeters')                                    %set paper size to figure size
set(currentFigure, 'PaperUnits', 'centimeters')
figureSize=get(currentFigure, 'Position');
set(currentFigure, 'PaperSize', figureSize(3:4))


%resultsTempSort=zeros(0,0)

temperatures=unique(results(2:end,2));
for i=1:length(temperatures);                                       %plot results for each temperature
    
    %resultsTempSort=[resultsTempSort cell2mat(results(strcmp(results(:,2),temperatures(i)),2+numStructures+1:2+2*numStructures))];
    
    currentFigure=figure;
    currentBar=bar3(cell2mat(results(strcmp(results(:,2),temperatures(i)),2+numStructures+1:2+2*numStructures)));
    for k = 1:length(currentBar)
        zdata = get(currentBar(k),'ZData');
        set(currentBar(k),'CData',zdata,...
            'FaceColor','interp')
    end
    set(gca,'CLimMode','auto')
    view(0, 90);
    title(temperatures(i));
    set(gca, 'XTickLabel',structureNames)
    set(gca, 'YTickLabel',num(strcmp(results(:,2),temperatures(i))))
    xlabel('structure')
    ylabel('parameter Set')
    colorbar
end