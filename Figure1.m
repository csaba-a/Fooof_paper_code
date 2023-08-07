%%%% Figure 1 %%%%%
clear all
close all
%% paths
define_paths
%%  features
type_power_spectrum='aperiodic'; %complete, periodic, aperiodic

% parcellation scheme
parc=4;

%define frequency bands
freq_bands = [1 4; %Delta
              4 8; %Theta
              8 13; %Alpha
              13 30]; %Beta

%% Load fooof output 
[normative_table,~]=load_fooof_data(fooof_output_path);

%% Choose the approprite ROIs based on parcellation scheme
if parc==2
    rois=normative_table.ROI_2;
elseif parc==4
    rois=normative_table.ROI_4;
end
%% Load and reorder Complete data
switch type_power_spectrum
    case 'complete'
        load_ieeg_psd %load ieeg complete data
        [MasterChannelTable]=reorder_ieeg(MasterChannelTable);
        remove_invalid_contacts % remove invalid contacts (e.g. in white matter etc)
        normative_table.complete_psd=MasterChannelTable.pxx_n(RAM_bool,:);

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

%% Calculate band powers
switch type_power_spectrum
    case 'complete' %Complete
        [rel_bp_complete,n_chan,n_bands]=calc_band_power(normative_table.complete_psd,freq_bands);

    case 'periodic'%Periodic
        [rel_bp_periodic,n_chan,n_bands]=calc_band_power(normative_table.flattened_psd,freq_bands);
end



%% Plot normative maps

switch type_power_spectrum
    case 'complete'
        %feature names
        feature_names={'Delta Bandpower','Theta Bandpower','Alpha Bandpower','Beta Bandpower'};
        % plot
        plot_normative_maps(rel_bp_complete, feature_names,rois,n_bands,colormap_idx,parc,analysis_location,figdir_1)
    case 'periodic'
        %feature names
        feature_names={'Periodic Delta Bandpower','Periodic Theta Bandpower','Periodic Alpha Bandpower','Periodic Beta Bandpower'};
        % plot
        plot_normative_maps(rel_bp_periodic, feature_names,rois,n_bands,colormap_idx,parc,analysis_location,figdir_1)
    case 'aperiodic'
        %feature names
        feature_names={'Aperiodic exponent'};
        n_bands=size(normative_table.aperiodic_cmps_2,2);
        % plot
        plot_normative_maps(normative_table.aperiodic_cmps_2, feature_names,rois,n_bands,colormap_idx_aperiodic,parc,analysis_location,figdir_1)

end