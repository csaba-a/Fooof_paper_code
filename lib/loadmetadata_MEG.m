function [metadata,MEG_auc]=loadmetadata_MEG(meg_data_path,MEG_auc)
%% LOADMETADATA_MEG
% Load the metadata and reorganize it to exclude TLE and pateints with
% insufficient imaging and resection info
% 
% Csaba Kozma
% CNNP Lab, Newcastle University
% July 2023
%%
MEG_auc = removevars(MEG_auc, "Var1");
metadata_full=readtable([meg_data_path,'meta_meg.csv']);
MEG_auc = sortrows(MEG_auc,"ID","ascend");
excluded_pats=metadata_full.ID(~ismember(metadata_full.ID,MEG_auc.ID));
metadata_full(ismember(metadata_full.ID,excluded_pats),:)=[];

ETLE_idx=find(strcmp(metadata_full.isNeocortical,'Yes'));
MEG_auc=MEG_auc(ETLE_idx,:);

metadata=table();
metadata.ID=MEG_auc.ID;
metadata.ILAE1=MEG_auc.outcome;
metadata.ILAE1(strcmp(metadata.ILAE1,'ILAE1'))={1};
metadata.ILAE1(strcmp(metadata.ILAE1,'ILAE2+'))={3};
metadata.ILAE1=cell2mat(metadata.ILAE1);
end