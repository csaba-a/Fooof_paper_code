function plot_normative_maps(normative_data,freq_bands, feature_names,rois,colormap_idx,parc,figdir,analysis_location,modality)
%% Calcualte mean values for normative maps and plot normative maps
% number of frequency bands
if strcmp(feature_names,'Aperiodic exponent')
    n_bands=size(normative_data,2);
else
    n_bands = size(freq_bands,1);
end
%% Load and reorder atlas
[atlas_tbl]=reorder_atlas();

% IF plotting MEG normative map, need to exclude subcortical areas
if strcmp(modality,'MEG')
    atlas_tbl.name{parc}(1:14)=[];
    atlas_tbl.xyz{parc}(1:14,:)=[];
end
%% Load brain surfaces and color scheme for plotting
load_brain_surface
normative_maps=nan(max(rois),n_bands);
for i=1:max(rois)
    locs=find(rois(:)==i);
    normative_maps(i,:)=mean(normative_data(locs,:),1);

end


%climits
climits=[nanmean(normative_maps,1)'-(nanstd(normative_maps,1))' nanmean(normative_maps,1)'+(nanstd(normative_maps,1))'];%[0.155 0.19];%Delta


for i=1:size(feature_names,2)
    bp_feat=feature_names{i};

    h=vis_norm_map_on_brain_for_fooof_paper(normative_maps(:,i),feature_names(i),atlas_tbl,'AtlasIndex',parc,'MarkerSize',200,'Colormap',sort(colormap_idx{i},'descend'),...
        'View',{'top','left','right'},'climits',climits(i,:),'PlotBrain',true,'FontSize',20);
    
    %Save
    %set(gcf,'renderer','painters');

    saveas(h, [figdir, bp_feat,'_parc_',num2str(parc),'_',modality,'.pdf']); % specify the full path
    close all
end

end