function [ errorSquared ] = generalErrorSquared2Fits( fitType1, xData1, yData1, parameterOrder1, fitType2, xData2, yData2, parameterOrder2, parameters )
%calculates squared error for two datasets x,y and two fit functions

%fittype1 first fit model 
%xdata1 independent variable data for first fit model
%ydata1 dependent variable data for first fit model
%parameterOrder1 order of values from parameters to supply to first fit

%fittype2 second fit model 
%xdata2 independent variable data for second fit model
%ydata2 dependent variable data for second fit model
%parameterOrder2 order of values from parameters to supply to second fit

%parameters all parameters for joint fit


n1=length(xData1);
n2=length(xData2);

mean1=mean(yData1);
mean2=mean(yData2);

inputParameters1=num2cell(parameters(parameterOrder1));
inputParameters2=num2cell(parameters(parameterOrder2));


errorSquared=0;


for i=1:n1
    errorSquared=errorSquared+(yData1(i)/(mean1)-fitType1(inputParameters1{:},xData1(i))/(mean1)).^2;
end

for i=1:n2
    errorSquared=errorSquared+(yData2(i)/(mean2)-fitType2(inputParameters2{:},xData2(i))/(mean2)).^2;
end


end

