%%%% Figure 1 %%%%%
clear all
close all
%% paths
define_paths
%%  features
type_power_spectrum='complete'; %complete, periodic, aperiodic

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
%% Load and reorder atlas
[atlas_tbl]=reorder_atlas();
%% Load brain surfaces and color scheme for plotting
load_brain_surface
%% colormap set up
% load color scheme
load('lib/colourmaps_yw.mat'); %this loads cmaps variable
colormap_idx=cell(5,1);
colormap_idx{1}=cmaps.delta_map; %Delta
colormap_idx{2}=cmaps.theta_map; %Theta
colormap_idx{3}=cmaps.alpha_map; %Alpha
colormap_idx{4}=cmaps.beta_map; %Beta
map = brewermap(256,'Oranges'); %Aperiodic (generate new color)
map(1,:) = [1 1 1];
colormap_idx{5}=map;

%% Calculate band powers
switch type_power_spectrum
    case 'complete' %Complete
        [rel_bp_complete,n_chan,n_bands]=calc_band_power(normative_table.complete_psd,freq_bands);

    case 'periodic'%Periodic
        [rel_bp_periodic,n_chan,n_bands]=calc_band_power(normative_table.flattened_psd,freq_bands);
end
%% Plotting preparation

% Take the mean across all controls to generate a whole brain map
prep_normative_map_data

%% Plot normative maps

plot_normative_maps