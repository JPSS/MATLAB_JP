%% select data
%selects data to open
%filename{N} are filenames of files
%pathname is pathnames of files

clearvars

[filename pathname]=uigetfile('*','feed me cluster data human','MultiSelect','on');


%% load data and analyze
%opens file, copies data into parameters and currentData, analyzes,
%then opens next file
%parameters are cells with parameter strings
%currentData is array of traceData x structures, with empty lines between
%traces omited

numStructures=15;                                                    %number of structures simulated
timeSteps=7200+1;                                                   %number of steps printed
traces=10;                                                           %number of traces per structure and parameter set

numParameterLines=27;                                               %number of lines after structure names before data starts
temperatureParameterLine=20;                                        %location of temperature value in beginning output
parametersText=cell(1+numStructures+numParameterLines,1);               %saves parameters text, each cell one line

numParameters=10;                                                   %number of parameters in filename
                

%structureNames={'42hbv3'};                                          %shorthand names of structures
structureNames={'v1','v2','v3','v4','v5','v6','v7','v8','v11','v13','v14','v15','ssp','RRv3','RR'};
%structureNames={'v1','v2','v3','v6','v7','v8','v13','v14','ssp','RR'};


numFiles=size(filename,2);                                          %number of files loaded
resultParameters=cell(numFiles+1,numParameters);                    %actual run parameters of each simulation file          

if ~(iscell(filename))                                              %check if a single file was loaded, if yes, adjust data type of filename
    numFiles=1;
    filename={filename};
end

currentData=zeros(timeSteps*traces,numStructures+1);                %temporary storage current file
results=cell(numFiles+1,1+1+2*numStructures);                         %analysis results of data


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
    
    %analysis part
                                                          
    results{i+1,1}=filename{i};                                     %put filenames in results
    
                                                                    %load parameters from file name
    [parameterStarts, parameterEnds]=regexp(results{i+1,1},'[0123456789.-]*');
    for j=1:length(parameterStarts)
        resultParameters{i+1,j}=results{i+1,1}(parameterStarts(j):parameterEnds(j));
    end
    
    tempLocation=strfind(filename{i}, '_t');                        %read temperature from filename
    results{i+1,2}=filename{i}(tempLocation+2:tempLocation+3);
    
    for j=1:numStructures                                           %set number of folded structures to max
        results{i+1,j+2}=traces;
        results{i+1,j+2+numStructures}=0;
    end
    
    for j=1:traces                                                  
        for k=1:numStructures

            if sum(currentData(j*timeSteps-30:j*timeSteps,k+1))~=0  %reduce number of folded structures if structure is unfolded at end of trace
                results{i+1,k+2}=results{i+1,k+2}-1;                %add constant non-fold time to average time of
                                                                    %folding
                results{i+1,k+2+numStructures}=results{i+1,k+2+numStructures}+currentData(timeSteps);
            else                                                    %find time of folding
                currentTimeStep=j*timeSteps;
                tempCurrentData=currentData((j-1)*timeSteps+1:j*timeSteps,[1,1+k]);
                tempCurrentData=tempCurrentData(tempCurrentData(:,2)~=0,1);
                                                                    %add time of folding to average time of folding
                results{i+1,k+2+numStructures}=results{i+1,k+2+numStructures}+tempCurrentData(end);
            end
            
        end
    end
    
    for j=1:numStructures                                           %calculate average time of folding
        results{i+1,j+2+numStructures}=results{i+1,j+2+numStructures}/traces;
    end
    
end

                                                                    %write column names from parameters text block
results{1,1}='filename';
results{1,2}='temperature';

for i=1:numStructures                                               %read names of structure files analyzed, write to results
    tempLocation=strfind(parametersText{i+1}, '/');
    results{1,i+2}=parametersText{i+1}(tempLocation(size(tempLocation,2))+1:size(parametersText{i+1},2));
    results{1,i+2+numStructures}=parametersText{i+1}(tempLocation(size(tempLocation,2))+1:size(parametersText{i+1},2));
end
    

%% save results to  results.txt

currentFile=fopen([pathname filesep 'results.txt'],'w');

for i=1:numFiles+1
    for j=1:2*numStructures+2
        if ischar(results{i,j})
            fprintf(currentFile,'%s\t',results{i,j});
        else
            fprintf(currentFile,'%f\t',results{i,j});
        end
    end
    fprintf(currentFile,'\n');
end

fclose(currentFile);


%% visualisation

for i=1:numParameters
    parameterTypes{i}=unique({resultParameters{2:end,i}})
    numParameterType{i}=length(parameterTypes{i})
end

                                                                %calculate avg fold time of each structure for one given parameter
currentFigure=figure;
currentAxes=gca;
for i=1:numParameters
    tempAvgs=zeros(numParameterType{i},numStructures);
    for j=1:numParameterType{i}
        for k=1:numStructures
            tempAvgs(j,k)=mean([results{strcmp(resultParameters(:,i),parameterTypes{i}(j)),2+numStructures+k}]);
        end
    end
    plot(tempAvgs)
    legend(structureNames)
    set(currentAxes, 'XTick',1:numParameterType{i}, 'XTickLabel',parameterTypes{i})
    pause
end
        




num=1:100;                                                       %parameter file counter

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