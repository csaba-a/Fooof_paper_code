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

%% Calculate group level AUC and pvalues
type_power_spectrum='complete'; %complete, periodic, aperiodic, max


% Load and reorder Complete data
switch type_power_spectrum
    case 'complete'
        % Plot AUC and beeswarm plots
        plot_group_auc(MEG_auc.drs_complete, metadata)
    case 'periodic'%Periodic
        % Plot AUC and beeswarm plots
        plot_group_auc(MEG_auc.drs_periodic, metadata)
    case 'aperiodic'%Aperiodic
        % Plot AUC and beeswarm plots
        plot_group_auc(MEG_auc.drs_aperiodic, metadata)
    case 'max'%max
        % Plot AUC and beeswarm plots
        plot_group_auc(MEG_auc.drs_max, metadata)
end

% Save plots
save_figs("Group_MEG",type_power_spectrum,analysis_location,figdir_4)

%% -- Supplementary -- %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compare abnormality distributions with correlation plots %%
%Run after generated all the abnormality matrices
% Separate good vs bad outcome pateints

good_outcome=MEG_auc(strcmp(MEG_auc.outcome,'ILAE1'),:);
bad_outcome=MEG_auc(strcmp(MEG_auc.outcome,'ILAE2+'),:);

% Define the type of DRS you want to compare to Complete Bp
% Drs
drs_to_compare = 'drs_max'; %drs_periodic, drs_aperiodic, drs_max

type = 'max'; %periodic, aperiodic, max


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
xlabel(['D_r_s score based on ', type, ' PSD'])
ylabel('D_r_s score based on Complete PSD')
ylim([0 1])
xlim([0 1])
set(gca, 'FontSize', 16)
title(['r=',num2str(round(correlat,2))])
set(gcf,'renderer','painters');
set(gca,'XTick',[])
set(gca,'YTick',[])
