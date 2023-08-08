%% %% Figure 4 %%%%%
clear all
close all

%% paths and features
define_paths

%% Load Drs values and metadata

%Drs values
[MEG_auc]=load_meg_aucs(meg_data_path);

%Metadata
[metadata,MEG_auc]=loadmetadata_MEG(meg_data_path,MEG_auc);

%% Plot group level AUC and pvalues
feature_labels=MEG_auc.Properties.VariableNames(contains(MEG_auc.Properties.VariableNames,'drs'));

for z=1:size(feature_labels,2)

    type=feature_labels{z}; %complete, periodic, aperiodic, max
    % Plot AUC and beeswarm plots
    plot_group_auc(MEG_auc.(type), metadata)
    % Save plots
    save_figs("Group_MEG",type,figdir_4)


end


%% -- Supplementary -- %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compare abnormality distributions with correlation plots %%
% Separate good vs bad outcome pateints
good_outcome=MEG_auc(strcmp(MEG_auc.outcome,'ILAE1'),:);
bad_outcome=MEG_auc(strcmp(MEG_auc.outcome,'ILAE2+'),:);


for z=1:size(feature_labels,2)
% Define the type of abnormality you want to compare to Complete Bp
% abnormality
%Skip complete as you always compare to this
if strcmp(feature_labels{z},'drs_complete')
    continue
end

drs_to_compare = feature_labels{z}; %select the auc labels

% Plotting
figure
set(gcf,'renderer','painters');
scatter(good_outcome.(drs_to_compare),good_outcome.drs_complete,100,'filled')
hold on
scatter(bad_outcome.(drs_to_compare),bad_outcome.drs_complete,100,'filled')

correlat=corr(MEG_auc.(drs_to_compare),MEG_auc.drs_complete,'rows','complete','Type','Spearman');
yline(0.5,'LineStyle','--')
xline(0.5,'LineStyle','--')
%txt = ['r=',num2str(round(correlat,2))];
%text(0.2,0.8,txt,'HorizontalAlignment','right','FontSize',16)
%legend({'Good Outcome','Bad Outcome'},'Location','northwest')
xlabel(['D_r_s score based on ', drs_to_compare(5:end), ' PSD'])
ylabel('D_r_s score based on complete PSD')
ylim([0 1])
xlim([0 1])
set(gca, 'FontSize', 16)
title(['r=',num2str(round(correlat,2))])
set(gcf,'renderer','painters');
set(gca,'XTick',[])
set(gca,'YTick',[])
%Save figure
saveas(gca, fullfile(figdir_4,strcat('MEG_Drs_scatter_Plot_',drs_to_compare,'_vs_drs_complete', '.pdf'))); % specify the full path
close all
end