% runtimes of simulations sorted by combination, runtime
endTimes;
startTimes=zeros(960,1);

runTimes=endTimes-startTimes;

%average runtime of ith combination
runTimesAvg=zeros(192,1);
runTimesStd=zeros(192,1);
for i=1:192
    runTimesAvg(i)=mean(runTimes(((i-1)*5+1):(i*5)));
    runTimesStd(i)=std(runTimes(((i-1)*5+1):(i*5)));
end


runTimes=endTimes-startTimes;

%original parameter sorting order for parameter combinations, first is cycled first
parameterMatrixOrg{1,1}='10000';
parameterMatrixOrg{2,1}='4';
parameterMatrixOrg{3,1}='56';
parameterMatrixOrg{3,2}='53';
parameterMatrixOrg{3,3}='50';
parameterMatrixOrg{3,4}='45';
parameterMatrixOrg{4,1}='-2.5340';
parameterMatrixOrg{4,2}='0.0000';
parameterMatrixOrg{5,1}='-8.3300';
parameterMatrixOrg{5,2}='-14.5222';
parameterMatrixOrg{6,1}='12.1722';
parameterMatrixOrg{6,2}='16.1592';
parameterMatrixOrg{7,1}='1.000';
parameterMatrixOrg{7,2}='0.003';
parameterMatrixOrg{7,3}='-0.003';
parameterMatrixOrg{8,1}='-0.002580';
parameterMatrixOrg{9,1}='-2.440';
parameterMatrixOrg{9,2}='-2.000';
parameterMatrixOrg{10,1}='5';

%combinationNumber count for each parameter value
combinationCountMatrix=ones(10,4);
parameterCounter=flipud(numParameterType');
for i=2:10
    parameterCounter{i}=parameterCounter{i}*parameterCounter{i-1};
end
parameterCounter
for i=2:10
    for j=1:numParameterType{11-i}
        combinationCountMatrix(i,j)=parameterCounter{i-1}*(j-1);
    end
end
combinationCountMatrix

%combination nr of ith loaded file
combinationNr=zeros(192,1);
for i=1:192
    for j=1:10
        combinationNr(i)=combinationNr(i)+combinationCountMatrix(strcmp(parameterMatrixOrg,resultParameters(i+1,j)));
    end
end

%average runtime of ith loaded file
runTimesAvgSorted=runTimesAvg(combinationNr);

plot(runTimesAvgSorted,[results{2:end,4}])
scatter(runTimesAvgSorted,[results{2:end,4}])

parameterSimulationTimes=zeros(10,4);
for i=1:10
    for j=1:numParameterType{11-i}
        parameterSimulationTimes(i,j)=mean(runTimesAvgSorted(strcmp(resultParameters(2:end,11-i),parameterMatrixOrg{i,j})));
    end
end

results((runTimesAvgSorted>1000000),1)
runTimesAvg((runTimesAvg>1000000))
runTimesStd((runTimesAvg>1000000))


temp=squeeze(averagedData);

surface(temp(2:30:7200,1:4:end))
surface(temp(2:30:7200,2:4:end))

surface(temp(2:3:720,1:4:end))
surface(temp(2:3:720,2:4:end))
surface(temp(2:3:720,3:4:end))
surface(temp(2:3:720,4:4:end))

surface(temp(2:3:720,1:4:24*4))
surface(temp(2:1:720,2:4:24*4),'EdgeColor','none','LineStyle','none')
surface(temp(2:3:720,3:4:24*4))
surface(temp(2:3:720,4:4:24*4))

[ results, averagedData ] = cluster_data_analysis(1,5,7201,{'v3'});

[ results, averagedData ] = cluster_data_analysis(2,2,7201,{'v1','v3'});

[ results, averagedData ] = cluster_data_analysis(10,10,7201);
[ results, averagedData ] = cluster_data_analysis(10,10,7201,{'v1','v2','v3','v6','v7','v8','v13','v14','ssp','RR'});

cluster_data_visualization( results, averagedData, 8, 11)
cluster_data_visualization( results, averagedData, 2, 11)

cluster_data_visualization( results, averagedData, 2, 8)

plot([results{5:7:end,14}])
hold on
plot([results{5:7:end,15}])
