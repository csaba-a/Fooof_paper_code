%% Run this script before all the others to define paths
analysis_location='/Users/c2056366/Documents/Fooof_paper_material/fooof_code_gms/Fooof_paper_code'; %Path to main library
fieldtrip_location='/Users/c2056366/Documents/fieldtrip-20221005'; %Fieldtrip
%define paths
addpath(genpath(analysis_location)) %custom scripts
cd(analysis_location);
% make figure path
cd ..
figure_location = pwd;

addpath(genpath(fieldtrip_location)) %Fieldtrip
%remove fieldtrip compat folder from path to plot
rmpath(genpath([fieldtrip_location,'/compat']));

complete_datapath='data'; %iEEG data
fooof_output_path='data/';%fooof output
metadata_path='data/metadata.mat'; %corresponding metadata
if ~exist([figure_location,'/Figures'])
    mkdir([figure_location,'/Figures'])
end
cd([figure_location,'/Figures'])
figdirs={'Figure1/','Figure2/','Figure3/','Figure4/'};

for i=1:length(figdirs)
if ~exist(figdirs{i})
    mkdir(figdirs{i})
end
end

figdir_1=[figure_location,'/Figures/Figure1/'];% figure save
figdir_2=[figure_location,'/Figures/Figure2/'];% figure save
figdir_3=[figure_location,'/Figures/Figure3/'];% figure save
figdir_4=[figure_location,'/Figures/Figure4/'];% figure save

%back to analysis location
cd(analysis_location); %back to main directory

meg_data_path = 'data/'; %meg data folder

