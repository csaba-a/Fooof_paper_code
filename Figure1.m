%%%% Figure 1 %%%%%
clear all
close all
%% paths
define_paths
%%  features
type_power_spectrum='aperiodic'; %complete, periodic, aperiodic

% parcellation scheme
parc=2;

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
        [MasterChannelTable,freq_needed]=load_ieeg_psd; %load ieeg complete data
        [MasterChannelTable]=reorder_ieeg(MasterChannelTable);
        [UCLH_bool,RAM_bool,MasterChannelTable]= remove_invalid_contacts(MasterChannelTable,freq_needed); % remove invalid contacts (e.g. in white matter etc)
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

%% Calculate band powers and define feature names
switch type_power_spectrum
    case 'complete' %Complete
        [normative_data,~]=calc_band_power(normative_table.complete_psd,freq_bands);
        %feature names
        feature_names={'Delta Bandpower','Theta Bandpower','Alpha Bandpower','Beta Bandpower'};

    case 'periodic'%Periodic
        [normative_data,~]=calc_band_power(normative_table.flattened_psd,freq_bands);
        %feature names
        feature_names={'Periodic Delta Bandpower','Periodic Theta Bandpower','Periodic Alpha Bandpower','Periodic Beta Bandpower'};

    case 'aperiodic' %Aperiodic exponent
        %feature names
        feature_names={'Aperiodic exponent'};
        normative_data=normative_table.aperiodic_cmps_2;
        colormap_idx=colormap_idx_aperiodic;

end



%% Plot normative maps


plot_normative_maps(normative_data,freq_bands, feature_names,rois,colormap_idx,parc,figdir_1,analysis_location)
