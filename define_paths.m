%% Run this script before all the others to define paths
analysis_location='/Users/c2056366/Documents/Fooof_paper_material/fooof_code_gms'; %Path to main library
fieldtrip_location='/Users/c2056366/Documents/fieldtrip-20221005'; %Fieldtrip

%define paths
addpath(genpath(analysis_location)) %custom scripts
cd(analysis_location);
addpath(genpath(fieldtrip_location)) %Fieldtrip
%remove fieldtrip compat folder from path to plot
rmpath(genpath([fieldtrip_location,'/compat']));

complete_datapath='data'; %iEEG data
fooof_output_path='data/';%fooof output
metadata_path='data/metadata.mat'; %corresponding metadata
if ~exist('Figures')
    mkdir('Figures')
end
cd('Figures')
figdirs={'Figure1/','Figure2/','Figure3/','Figure4/'};

for i=1:length(figdirs)
if ~exist(figdirs{i})
    mkdir(figdirs{i})
end
end

cd(analysis_location); %back to main directory
figdir_1='/Figures/Figure1/';% figure save
figdir_2='/Figures/Figure2/';% figure save
figdir_3='/Figures/Figure3/';% figure save
figdir_4='/Figures/Figure4/';% figure save

meg_data_path = 'data/'; %meg data folder

