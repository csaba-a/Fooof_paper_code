function [MEG_normative_table]=load_MEG_material(meg_data_path,parc)

%% Load BP values and misc
% Load original data
load([meg_data_path,'cardiff_bp_parc_',num2str(parc),'.mat']);

% Reshape original PSD array Cardiff
psd_collapsed = permute(cardiff_bp,[1 3 2]);
psd_collapsed = reshape(psd_collapsed,[],size(cardiff_bp,2),1);
cardiff_bp_original = rmmissing(psd_collapsed);

clear Cardiff_bp UCL_bp psd_collapsed cardiff_bp

% Load Periodic data
load([meg_data_path,'cardiff_bp_flattened_30_parc_',num2str(parc),'.mat']);

% Reshape periodic PSD array Cardiff
psd_collapsed = permute(cardiff_bp,[1 3 2]);
psd_collapsed = reshape(psd_collapsed,[],size(cardiff_bp,2),1);
cardiff_bp_periodic = rmmissing(psd_collapsed);

clear Cardiff_bp UCL_bp psd_collapsed cardiff_bp

% Load Aperiodic data
MEG_normative_table=readtable([meg_data_path,'Cardiff_aperiodic_30_parc_',num2str(parc),'.csv']);


% Organize all data into the Cradiff table
MEG_normative_table.complete_psd=cardiff_bp_original;
MEG_normative_table.periodic_psd=cardiff_bp_periodic;
% Load atlas for ROIs and xyz


end