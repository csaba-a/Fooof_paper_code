%%%% Figure 1 %%%%%
clear all
close all
%% paths
%% Load paths
analysis_location='/Users/c2056366/Documents/Fooof_paper_material/fooof_code_gms/Fooof_paper_code'; %Path to main library
figdir='/Users/c2056366/Documents/Fooof_paper_material/fooof_code_gms/Figures/Supplementary/';% figure save
meg_data_path = '/Users/c2056366/Documents/Fooof_paper_material/fooof_code_gms/Fooof_paper_code/data/Supplementary/';
%%  features
type_power_spectrum='aperiodic'; %complete, periodic, aperiodic
modality='MEG';
% parcellation scheme
parc=2;

%define frequency bands
freq_bands = [1 4; %Delta
              4 8; %Theta
              8 13; %Alpha
              13 30]; %Beta

%% Load MEG data and atlas and surfaces
[MEG_normative_table]=load_MEG_material(meg_data_path,parc);

%% Choose the approprite ROIs based on parcellation scheme
if parc==2
    rois=MEG_normative_table.ROI_2;
elseif parc==4
    rois=MEG_normative_table.ROI_4;
end

%% colormap set up
% load color scheme
load('lib/colourmaps_yw.mat'); %this loads cmaps variable
colormap_idx=cell(4,1);
colormap_idx{1}=cmaps.delta_map; %Delta
colormap_idx{2}=cmaps.theta_map; %Theta
colormap_idx{3}=cmaps.alpha_map; %Alpha
colormap_idx{4}=cmaps.beta_map; %Beta
map = brewermap(256,'Oranges'); %Aperiodic (generate new color)
map(1,:) = [1 1 1];
colormap_idx_aperiodic=mat2cell(map,size(map,1),size(map,2));

%% Calculate band powers and define feature names
switch type_power_spectrum
    case 'complete' %Complete
        normative_data=MEG_normative_table.complete_psd;
        %feature names
        feature_names={'Delta Bandpower','Theta Bandpower','Alpha Bandpower','Beta Bandpower'};

    case 'periodic'%Periodic
       normative_data=MEG_normative_table.periodic_psd;
        %feature names
        feature_names={'Periodic Delta Bandpower','Periodic Theta Bandpower','Periodic Alpha Bandpower','Periodic Beta Bandpower'};
        
    case 'aperiodic' %Aperiodic exponent
        %feature names
        feature_names={'Aperiodic exponent'};
        normative_data=MEG_normative_table.aperiodic_cmps_2;
        colormap_idx=colormap_idx_aperiodic;

end


%% Plot normative maps


plot_normative_maps(normative_data,freq_bands, feature_names,rois,colormap_idx,parc,figdir,analysis_location, modality)
