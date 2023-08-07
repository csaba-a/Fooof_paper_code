
function [resected,spared]=plot_abnormalities(UCLROIdata_abnormality,UCLROIcontainsresected,resectedThresh,parc)
%% plot max z score
% Plotting abnormality for individual patients with brain surface from top
% view
% Csaba Kozma
% CNNP Lab, Newcastle University
% July 2023
resected=find(UCLROIcontainsresected>resectedThresh);
spared=find(UCLROIcontainsresected<=resectedThresh & UCLROIcontainsresected~=-1);

%% Load and reorder atlas
[atlas_tbl]=reorder_atlas();

%% load colormaps
load('lib/colourmaps_yw.mat'); %this loads cmaps variable

%% Plotting
vis_norm_map_on_brain_for_fooof_paper(UCLROIdata_abnormality,{'Abnormality'},atlas_tbl,'AtlasIndex',parc,'MarkerSize',400,'Colormap',sort(cmaps.z_score_map,'descend'),...
                'View',{'top'},'climits',[0 4],'PlotBrain',true,'FontSize',20,'resected',resected);

end
