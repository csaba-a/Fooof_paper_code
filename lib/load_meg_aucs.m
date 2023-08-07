function [MEG_auc]=load_meg_aucs(auc_data_path)
%% LOAD_MEG_AUCS
% load the already calculated Drs scores for every patients across complete
% band power, periodic band power, aperiodic exponent and max abnormality
% approaches
% 
% Csaba Kozma
% CNNP Lab, Newcastle University
% July 2023
%%
MEG_auc=readtable([auc_data_path,'MEG_auc_cs_complete_30.csv']);
MEG_auc_flattaned=readtable([auc_data_path,'MEG_auc_cs_flat_30.csv']);
MEG_auc_aperiodic=readtable([auc_data_path,'MEG_auc_cs_slopes_30.csv']);
MEG_auc_max=readtable([auc_data_path,'MEG_auc_cs_mixed.csv']);

MEG_auc.drs_periodic=MEG_auc_flattaned.drs;
MEG_auc.drs_aperiodic=MEG_auc_aperiodic.drs;
MEG_auc.drs_max=MEG_auc_max.drs;
MEG_auc.Properties.VariableNames(4) = "drs_complete";
end