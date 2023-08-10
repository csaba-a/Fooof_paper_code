function [metadata,nbUCLsubj]=load_metadata(metadata_path)
%% LOAD METADATA: loads metadata for the iEEG cohort and selects out the ones needed for analysis
% Csaba Kozma
% CNNP Lab, Newcastle University
% July 2023
%% Set up the metadata including outcome and ID
IDP={};ILAE1=[];
m=table(IDP,ILAE1);
clear metadata
load(metadata_path)
startlocs=size(m.IDP,1)+1;
endlocs=startlocs+(size(metadata,1)-1);
m.IDP(startlocs:endlocs)=convertStringsToChars(metadata.IDP);
m.ILAE1(startlocs:endlocs)=metadata.ILAE1;


metadata=m;
clear m startlocs endlocs
nbUCLsubj=size(metadata,1);
end