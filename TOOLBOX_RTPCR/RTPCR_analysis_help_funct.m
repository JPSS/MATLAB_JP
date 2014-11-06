function [ results ] = RTPCR_analysis_help_funct( reactionPartnersSpeciesId, data, fits )
%fits rt pcr data using model:
%probabilitiy of A+B reaction is pA*pB
%returns distance vector between data vector and fit vector

results=data-fits(reactionPartnersSpeciesId(1,:)).*fits(reactionPartnersSpeciesId(2,:));

end

