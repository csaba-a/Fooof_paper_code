%%%% Figure 3 %%%%%
clear all
close all

%% paths and features
define_paths

%% Run this section once
% Load and set up metadata and select example patient
[metadata,~]=load_metadata(metadata_path);

% parcellation scheme
parc=2;
% Resection threshold
resectedThresh=0.25;
%define frequency bands
freq_bands = [1 4; %Delta
    4 8; %Theta
    8 13; %Alpha
    13 30]; %Beta

% Load fooof output and set up path and periodic/original data
[normative_table,patient_table]=load_fooof_data(fooof_output_path);


% Generate table for abnormalities,aucs and add outcome column
abnormalities_tbl=table();
auc_tbl=table();
auc_tbl.outcome = metadata.ILAE1; %ILAE 1 year outcome add to auc table


%% Repeate for all types -- Max can only be calculated if Periodic and Aperiodic are calculated


% Load and reorder Complete data

[MasterChannelTable,freq_needed]=load_ieeg_psd; %load ieeg complete data
[MasterChannelTable]=reorder_ieeg(MasterChannelTable);
[UCLH_bool,RAM_bool,MasterChannelTable]= remove_invalid_contacts(MasterChannelTable,freq_needed); % remove invalid contacts (e.g. in white matter etc)



%% Calculate band powers and abnromality for complete and periodic
%Calculate relative band power for Complete
[rel_bp_complete,n_chan]=calc_band_power(MasterChannelTable.pxx_n,freq_bands);
normative_table.rel_bp=rel_bp_complete(RAM_bool,:);
patient_table.rel_bp=rel_bp_complete(UCLH_bool,:);

%Calculate abnormality for complete bp
[abnormalities_tbl.complete_abn,~]=calc_abnormality(normative_table, patient_table,metadata.IDP,freq_bands,resectedThresh, parc);


%Calculate relative band power for Periodic
flattened_pxx=[normative_table.flattened_psd; patient_table.flattened_psd];
[rel_bp_periodic,~]=calc_band_power(flattened_pxx,freq_bands);
%Data source indices
data_source=[normative_table.DataSource; patient_table.DataSource];


%Define boolian for separating datasets
UCLH_bool = data_source=="UCLH";
RAM_bool = data_source=="RAM";
normative_table.rel_bp=rel_bp_periodic(RAM_bool,:);
patient_table.rel_bp=rel_bp_periodic(UCLH_bool,:);

%Calculate abnormality for periodic bp
[abnormalities_tbl.periodic_abn,~]=calc_abnormality(normative_table, patient_table,metadata.IDP,freq_bands,resectedThresh, parc);


% Calculate abnormality, plot and save abnormality plots

%% Calculate abnormality for Aperiodic exponent and Max across periodic bp and aperiodic exponent

%Aperiodic abnormality
[abnormalities_tbl.aperiodic_abn,UCLROIcontainsresected]=calc_abnormality_aperiodic(normative_table, patient_table,metadata.IDP,resectedThresh, parc);

%Max abnormality across periodic and aperiodic 
for z=1:length(metadata.IDP)
    abnormalities_tbl.max_abn(:,z)=max(abs(abnormalities_tbl.periodic_abn(:,z)),abs(abnormalities_tbl.aperiodic_abn(:,z)));
end
%% Calculate AUCs

for i=1:size(abnormalities_tbl,2)
    
    %Calculate AUC
    [auc_tbl.(abnormalities_tbl.Properties.VariableNames{i})]=calc_auc(abnormalities_tbl.(abnormalities_tbl.Properties.VariableNames{i}),UCLROIcontainsresected,metadata.IDP,resectedThresh);

    % Plot AUC and beeswarm plots
    plot_group_auc(auc_tbl.(abnormalities_tbl.Properties.VariableNames{i}), metadata)
    % Save plots
    save_figs("Group_AUC_",abnormalities_tbl.Properties.VariableNames{i},figdir_3)
end


%% %% -- Supplementary -- %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare Drs distributions with correlation plots %%
%Run after generated all the abnormality matrices
% Separate good vs bad outcome pateints

good_outcome=auc_tbl(auc_tbl.outcome<3,:);
bad_outcome=auc_tbl(auc_tbl.outcome>2,:);

% define feature labels
feature_labels=auc_tbl.Properties.VariableNames(contains(auc_tbl.Properties.VariableNames,'abn'));

for z=1:size(feature_labels,2)
% Define the type of abnormality you want to compare to Complete Bp
% abnormality
%Skip complete as you always compare to this
if strcmp(feature_labels{z},'complete_abn')
    continue
end

drs_to_compare = feature_labels{z}; %select the auc labels


% Plotting
figure
set(gcf,'renderer','painters');
scatter(good_outcome.(drs_to_compare),good_outcome.complete_abn,100,'filled')
hold on
scatter(bad_outcome.(drs_to_compare),bad_outcome.complete_abn,100,'filled')

correlat=corr(auc_tbl.(drs_to_compare),auc_tbl.complete_abn,'rows','complete','Type','Spearman');
yline(0.5,'LineStyle','--')
xline(0.5,'LineStyle','--')
%txt = ['r=',num2str(round(correlat,2))];
%text(0.2,0.8,txt,'HorizontalAlignment','right','FontSize',16)
%legend({'Good Outcome','Bad Outcome'},'Location','northwest')
xlabel(['D_r_s score based on ', drs_to_compare(1:end-4), ' PSD'])
ylabel('D_r_s score based on Complete PSD')
axis equal
ylim([0 1])
xlim([0 1])
set(gca, 'FontSize', 16)
title(['r=',num2str(round(correlat,2))])
set(gcf,'renderer','painters');
set(gca,'XTick',[])
set(gca,'YTick',[])
%Save figure
saveas(gca, fullfile(figdir_3,strcat('iEEG_Drs_scatter_Plot_',drs_to_compare,'_vs_drs_complete', '.pdf'))); % specify the full path
close all

end

