%% load RT PCR DATA pre-prepared
%each column is one reaction
%first lane: interaction partner 1
%second lane: interaction partner 2
%further lanes: "probability of bound state"

[filename pathname]=uigetfile('*.txt','THIS BETTER BE THE RIGHT FILE.txt');
cd(pathname)


%% read data

RTPCRData=dlmread([ pathname filename ]);

cd(matlab_dir);

numberOfTraces=size(RTPCRData,2)                        %number of traces loaded
numberOfPoints=size(RTPCRData,1)-2                      %number of datapoints per trace

%% set up equation system

reactionPartners=RTPCRData(1:2,1:numberOfTraces)        %reaction partners in each reaction
reactionPartnersSpeciesId=zeros(2,numberOfTraces);      %species ids of reaction partners in each reaction
species=unique(reactionPartners);                       %array of unique species
numberOfSpecies=size(species,1)                         %number of unique species

lowerBounds=zeros(1,numberOfSpecies);                   %lower bounds species probability
upperBounds=ones(1,numberOfSpecies);                    %upper bounds species probability

speciesProbability=ones(numberOfPoints,numberOfSpecies);  %resulting fits of species probabilities

for i=1:numberOfTraces                                  %species id partners in each reaction
    reactionPartnersSpeciesId(1,i)=find(species==reactionPartners(1,i));
    reactionPartnersSpeciesId(2,i)=find(species==reactionPartners(2,i));
end

for i=1:numberOfPoints                                  %fit species probabilities for each timestep
    if i>1                                              %set starting values to last fit results
        speciesProbability(i,:)=speciesProbability(i-1,:);
    end
                                                       
    speciesProbability(i,:)=lsqnonlin(@(fits)RTPCR_analysis_help_funct(reactionPartnersSpeciesId,RTPCRData(i+2,:),fits),speciesProbability(i,:),lowerBounds,upperBounds);
    
end


newFileName=inputdlg('filename species probability', 'filename species probability', 1, {'species probability.txt'});
dlmwrite([pathname newFileName{1}], speciesProbability, 'delimiter' , '\t');

newFileName=inputdlg('filename reaction probabilities', 'filename reaction probabilities', 1, {'reaction probabilities.txt'});
dlmwrite([pathname newFileName{1}], reactionProbabilities, 'delimiter' , '\t');    
        



