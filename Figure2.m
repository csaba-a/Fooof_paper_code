%%%% Figure 2 %%%%%
clear all
close all

%% paths and features
define_paths

%% Rerun this section when assign new example patient
% Example patient
example_patient="1216";%1216 --good outcome (Patient1) ; 910 -- bad outcome (Patient2)


% Load and set up metadata and select example patient
[metadata,~]=load_metadata(metadata_path);

%select out the right example patient to analyse and plot
metadata(metadata.IDP~=example_patient,:)=[];

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


% Generate table for abnormalities
abnormalities_tbl=table();

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



%% Calculate abnormality for Aperiodic exponent and Max across periodic bp and aperiodic exponent

%Aperiodic abnormality
[abnormalities_tbl.aperiodic_abn,UCLROIcontainsresected]=calc_abnormality_aperiodic(normative_table, patient_table,metadata.IDP,resectedThresh, parc);

%Max abnormality across periodic and aperiodic 
for z=1:length(metadata.IDP)
    abnormalities_tbl.max_abn(:,z)=max(abs(abnormalities_tbl.periodic_abn(:,z)),abs(abnormalities_tbl.aperiodic_abn(:,z)));
end



%% Plot abnormality brain plots
for i=1:size(abnormalities_tbl,2)

    [resected,spared]=plot_abnormalities(abnormalities_tbl.(abnormalities_tbl.Properties.VariableNames{i}),UCLROIcontainsresected,resectedThresh,parc);
    % Plot DRS plots
    drs_plot_indiv(resected,spared,abnormalities_tbl.(abnormalities_tbl.Properties.VariableNames{i}));
    % Save plots
    save_figs(example_patient,abnormalities_tbl.Properties.VariableNames{i},figdir_2)
end
%% Compare abnormality distributions with correlation plots %%

for z=1:size(abnormalities_tbl,2)
% Define the type of abnormality you want to compare to Complete Bp
% abnormality
%Skip complete as you always compare to this
if strcmp(abnormalities_tbl.Properties.VariableNames{z},'complete_abn')
    continue
end

abnormality_to_compare = abnormalities_tbl.Properties.VariableNames{z}; %abnormalities_tbl.aperiodic_abn or abnormalities_tbl.max_abn

% Plotting

figure


set(gcf,'renderer','painters');
scatter(abnormalities_tbl.(abnormality_to_compare)(spared),abnormalities_tbl.complete_abn(spared),200,'filled')
hold on
scatter(abnormalities_tbl.(abnormality_to_compare)(resected),abnormalities_tbl.complete_abn(resected),200,'filled')

correlat=corr(abnormalities_tbl.(abnormality_to_compare),abnormalities_tbl.complete_abn,'rows','complete','Type','Spearman');
%legend({'Spared','Resected'},'Location','northwest')
xlabel([abnormality_to_compare(1:end-4),' abnormality'])
ylabel('complete abormality')
set(gca, 'FontSize', 16)
title(['r=',num2str(round(correlat,2))])
axis equal
xlim([0 4])
ylim([0 4])
set(gca,'XTick',[0,2,4])
set(gca,'YTick',[0,2,4])

%Save figure
saveas(gca, fullfile(figdir_2,strcat(example_patient,'_scatter_Plot_',abnormality_to_compare,'_vs_complete_abnormality', '.pdf'))); % specify the full path
close all
end