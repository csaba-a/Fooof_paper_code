%%%% Supp. MEG normative maps%%%%

clear all
close all
%% Load paths
addpath(genpath('/Users/c2056366/Documents/Fooof_paper_material')) %custom scripts
addpath(genpath('/Users/c2056366/Documents/fieldtrip-20221005')) %Fieldtrip
figdir='/Users/c2056366/Documents/Fooof_paper_material/Figures/Supplementary/';% figure save
meg_data_path = '/Users/c2056366/Documents/Fooof_paper_material/fooof_output/';

%% Complete Band powe or Periodic Band power
type_power_spectrum='periodic';
parc=2;
%% colormap set up
% load color scheme
load('lib/colourmaps_yw.mat');

colormap_idx=cell(5,1);
colormap_idx{1}=cmaps.delta_map;
colormap_idx{2}=cmaps.theta_map;
colormap_idx{3}=cmaps.alpha_map;
colormap_idx{4}=cmaps.beta_map;
map = brewermap(256,'Oranges'); 
map(1,:) = [1 1 1];
colormap_idx{5}=map;
%% Load MEG data and atlas and surfaces
[atlas,cardiff_bp_original,cardiff_bp_periodic,cardiff_table]=load_MEG_material(meg_data_path,parc);

%% Choose the approprite ROIs based on parcellation scheme
if parc==2
    rois=cardiff_table.ROI_2;
elseif parc==4
    rois=cardiff_table.ROI_4;
end

%% Plotting prep organize data

% Take the mean across all controls to generate a whole brain map

switch type_power_spectrum
    case 'complete'

        %Complete PSD
        
        normative_maps_BP_psd=cell(5,1);
        for i=1:max(rois)
            locs=find(rois(:)==i);
            normative_maps_BP_psd{1}(i)=mean(cardiff_bp_original(locs,1));
            normative_maps_BP_psd{2}(i)=mean(cardiff_bp_original(locs,2));
            normative_maps_BP_psd{3}(i)=mean(cardiff_bp_original(locs,3));
            normative_maps_BP_psd{4}(i)=mean(cardiff_bp_original(locs,4));
            normative_maps_BP_psd{5}(i)=mean(cardiff_table.aperiodic_cmps_2(locs));

        end

    case 'periodic'
        % flat PSD

        normative_maps_BP_flat_psd=cell(4,1);
        for i=1:max(rois)
            locs=find(rois(:)==i);
            normative_maps_BP_flat_psd{1}(i)=mean(cardiff_bp_periodic(locs,1));
            normative_maps_BP_flat_psd{2}(i)=mean(cardiff_bp_periodic(locs,2));
            normative_maps_BP_flat_psd{3}(i)=mean(cardiff_bp_periodic(locs,3));
            normative_maps_BP_flat_psd{4}(i)=mean(cardiff_bp_periodic(locs,4));


        end

end

%% PLOT

switch type_power_spectrum
    % for Complete PSD

    case 'complete'

        %climits
        climits=nan(5,2);
        climits(1,:)=[0.12 0.30]; %Delta
        climits(2,:)=[0.13 0.220]; %Theta
        climits(3,:)=[0.21 0.38];%Alpha
        climits(4,:)=[0.20 0.30];%Beta
        climits(5,:)=[0.4 0.9]; %Aperiodic exponent

        %feature names
        feat_names={'Delta Bandpower','Theta Bandpower','Alpha Bandpower','Beta Bandpower','Aperiodic'};

        for i=1:size(feat_names,2)
            bp_feat=feat_names{i};

            h=vis_norm_map_on_brain_for_fooof_paper(normative_maps_BP_psd{i},feat_names(i),atlas,'AtlasIndex',parc,'MarkerSize',200,'Colormap',sort(colormap_idx{i},'descend'),...
                'View',{'left','right'},'climits',climits(i,:),'PlotBrain',true,'FontSize',20);

            %Save
            %set(gcf,'renderer','painters');

            saveas(h, [figdir, bp_feat,'MEG.pdf']); % specify the full path
            close all
        end

    case 'periodic'
        % PLOT flat PSD
        %climits
        climits=nan(5,2);
        climits(1,:)=[0.05 0.1];%[min(normative_maps_BP_flat_psd{i}) max(normative_maps_BP_flat_psd{i})]; %Delta
        climits(2,:)=[0.1 0.17];%[min(normative_maps_BP_flat_psd{i}) max(normative_maps_BP_flat_psd{i})]; %Theta
        climits(3,:)=[0.23 0.35];%[min(normative_maps_BP_flat_psd{i}) max(normative_maps_BP_flat_psd{i})] %Alpha
        climits(4,:)=[0.40 0.50];%[min(normative_maps_BP_flat_psd{i}) max(normative_maps_BP_flat_psd{i})]; %Beta

        %feature names
        feat_names={'Periodic Delta Bandpower','Periodic Theta Bandpower','Periodic Alpha Bandpower','Periodic Beta Bandpower'};

        for i=1:size(feat_names,2)
            bp_feat=feat_names{i};

             h=vis_norm_map_on_brain_for_fooof_paper(normative_maps_BP_flat_psd{i},feat_names(i),atlas,'AtlasIndex',parc,'MarkerSize',200,'Colormap',sort(colormap_idx{i},'descend'),...
                'View',{'left','right'},'climits',climits(i,:),'PlotBrain',true,'FontSize',20);

            %Save
            set(gcf,'renderer','painters');
            saveas(h, [figdir,bp_feat, 'MEG_periodic.pdf']); % specify the full path
            close all
        end
end