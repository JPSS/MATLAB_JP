function [ ] = cluster_data_visualization( results, averagedData, param1, param2)
%% plot analyzed cluster data

numParameters=10;                                                   %number of parameters in filename
numDataPoints=size(averagedData,1);
numStructures=size(averagedData,2);
numFiles=size(results,1)-1;

structureNames=results(1,2+numParameters:1+numParameters+numStructures);

resultsIndividual=cell(numFiles*numStructures,numParameters+2);
for i=0:numFiles-1
    for j=1:numStructures
        resultsIndividual(i*numStructures+j,1:end-2)=results(i+2,2:1+numParameters);
        resultsIndividual(i*numStructures+j,end-1)=structureNames(j);
        resultsIndividual(i*numStructures+j,end)=results(i+2,1+numParameters+numStructures+j);
    end
end

parameterTypes=cell(1,numParameters+1);
numParameterType=cell(1,numParameters+1);
for i=1:numParameters+1
    parameterTypes{i}=unique(resultsIndividual(:,i));
    numParameterType{i}=length(parameterTypes{i});
end

parameterFoldingTime=zeros(numParameterType{param1},numParameterType{param2});

param3=8;

for i=1:numParameterType{param1}

    locationsParam1=strcmp(resultsIndividual(:,param1),parameterTypes{param1}(i));
    for j=1:numParameterType{param2}

        locationsParam2=strcmp(resultsIndividual(:,param2),parameterTypes{param2}(j));
        
        for k=1:numParameterType{param3}
            locationsParam3=strcmp(resultsIndividual(:,param3),parameterTypes{param3}(k));
        
            parameterFoldingTime(i,j,k)=mean([resultsIndividual{locationsParam1 & locationsParam2 & locationsParam3,2+numParameters}]);
        end
    
    end
end

currentFigure=figure;       
for i=1:numParameterType{param3}
parameterTypes{param3}(i)
currentAxes=gca;
plot(parameterFoldingTime(:,:,i))

set(0,'defaultAxesLineStyleOrder','-|--|:')
ylabel('folding time [sec]')
axis([-inf inf 0 numDataPoints])
set(currentAxes, 'XTick',1:numParameterType{param1}, 'XTickLabel',parameterTypes{param1})
legend(parameterTypes{param2})
set(currentAxes, 'FontSize',20)
pause
end

return

currentFigure=figure;       
for i=1:numStructures
    clf
    surface(squeeze(averagedData(2:300,i,:)),'FaceColor','interp');
    currentAxes=gca;
    title(results{1,1+numParameters+i});
    set(currentAxes, 'XTick',1:numParameterType{param1}, 'XTickLabel',results(2:end,1+param1))
    pause
end

for i=1:numParameterType{param1}
    clf
    surface(squeeze(averagedData(2:300,:,i)),'FaceColor','interp');
    currentAxes=gca;
    title(results{1+i,1+param1});
    set(currentAxes, 'XTick',1:numStructures, 'XTickLabel',{results{1,2+numParameters:1+numParameters+numStructures}})
    pause
end

%surface(temp(2:1:720,2:4:24*4),'EdgeColor','none','LineStyle','none')


return
                                                                %calculate avg fold time of each structure for one given parameter
temperatures=unique(results(2:end,1+temperatureNameLocation));
                                                                
currentFigure=figure;       
currentAxes=gca;
for i=1:numParameters
    tempAverages=zeros(numParameterType{i},length(temperatures));
    for j=1:numParameterType{i}
        for k=1:length(temperatures)
            tempAverages(j,k)=mean([results{strcmp(results(2:end,i+1),parameterTypes{i}(j)) & strcmp(results(2:end,1+temperatureNameLocation),temperatures(k)),4}]);
        end
    end
    plot(transpose(tempAverages),'LineWidth',5) %plot fold time over temperature for one constant parameter averaged over all other parameters
    legend(parameterTypes{i})
    ylabel('folding time [sec]')
    xlabel('Temperature [C]')
    axis([-inf inf 0 size(averagedData,1)])
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




return

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



return