%%%% Figure 3 %%%%%
clear all
close all

%% paths and features
define_paths

%% Run this section once
% Load and set up metadata and select example patient
[metadata,UCLsubjID,~]=load_metadata(metadata_path);

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


% Generate table for abnormalities and add outcome column
auc_tbl=table();

auc_tbl.outcome = metadata.ILAE1; %ILAE 1 year outcome add to auc table


%% Repeate for all types -- Max can only be calculated if Periodic and Aperiodic are calculated

type_power_spectrum='max'; %complete, periodic, aperiodic, max

%

% Load and reorder Complete data
switch type_power_spectrum
    case 'complete'
        load_ieeg_psd %load ieeg complete data
        [MasterChannelTable]=reorder_ieeg(MasterChannelTable);
        remove_invalid_contacts % remove invalid contacts (e.g. in white matter etc)

end

% Calculate band powers
switch type_power_spectrum
    case 'complete' %Complete
        [rel_bp,n_chan,n_bands]=calc_band_power(MasterChannelTable.pxx_n,freq_bands);
        normative_table.rel_bp=rel_bp(RAM_bool,:);
        patient_table.rel_bp=rel_bp(UCLH_bool,:);
    case 'periodic'%Periodic
        %Combine data for relative bp calculation
        flattened_pxx=[normative_table.flattened_psd; patient_table.flattened_psd];
        %Data source indices
        data_source=[normative_table.DataSource; patient_table.DataSource];
        [rel_bp,n_chan,n_bands]=calc_band_power(flattened_pxx,freq_bands);

        %Define boolian for separating datasets
        UCLH_bool = data_source=="UCLH";
        RAM_bool = data_source=="RAM";
        normative_table.rel_bp=rel_bp(RAM_bool,:);
        patient_table.rel_bp=rel_bp(UCLH_bool,:);

end

% Calculate abnormality, plot and save abnormality plots

switch type_power_spectrum
    case 'complete' %Complete
        [UCLROIdata_abnormality_complete,UCLROIcontainsresected]=calc_abnormality(normative_table, patient_table,UCLsubjID,n_bands,resectedThresh, parc);
        [aucdrs_complete]=calc_auc(UCLROIdata_abnormality_complete,UCLROIcontainsresected,UCLsubjID,resectedThresh);
        auc_tbl.complete_auc=aucdrs_complete;
        % Plot AUC and beeswarm plots
        plot_group_auc(aucdrs_complete, metadata)

    case 'periodic' %Periodic
        [UCLROIdata_abnormality_periodic,UCLROIcontainsresected]=calc_abnormality(normative_table, patient_table,UCLsubjID,n_bands,resectedThresh, parc);
        [aucdrs_periodic]=calc_auc(UCLROIdata_abnormality_periodic,UCLROIcontainsresected,UCLsubjID,resectedThresh);
        auc_tbl.periodic_auc=aucdrs_periodic;
        % Plot AUC and beeswarm plots
        plot_group_auc(aucdrs_periodic, metadata)

    case 'aperiodic' %Aperiodic
        [UCLROIdata_abnormality_aperiodic,UCLROIcontainsresected]=calc_abnormality_aperiodic(normative_table, patient_table,UCLsubjID,resectedThresh, parc);
        [aucdrs_aperiodic]=calc_auc(UCLROIdata_abnormality_aperiodic,UCLROIcontainsresected,UCLsubjID,resectedThresh);
        auc_tbl.aperiodic_auc=aucdrs_aperiodic;
        % Plot AUC and beeswarm plots
        plot_group_auc(aucdrs_aperiodic, metadata)
    case 'max' %Max abnormality across periodic and aperiodic
        for i=1:size(UCLROIdata_abnormality_periodic,1)
            for z=1:length(UCLsubjID)
                UCLROIdata_abnormality_max(i,z)=max(abs(UCLROIdata_abnormality_periodic(i,z)),abs(UCLROIdata_abnormality_aperiodic(i,z)));
            end
        end
        [aucdrs_max]=calc_auc(UCLROIdata_abnormality_max,UCLROIcontainsresected,UCLsubjID,resectedThresh);
        auc_tbl.max_auc=aucdrs_max;
        % Plot AUC and beeswarm plots
        plot_group_auc(aucdrs_max, metadata)
end

% Save plots
save_figs("Group",type_power_spectrum,analysis_location,figdir_3)



%% %% -- Supplementary -- %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare abnormality distributions with correlation plots %%
%Run after generated all the abnormality matrices
% Separate good vs bad outcome pateints

good_outcome=auc_tbl(auc_tbl.outcome<3,:);
bad_outcome=auc_tbl(auc_tbl.outcome>2,:);

% Define the type of DRS you want to compare to Complete Bp
% Drs
drs_to_compare = 'max_auc'; %periodic_auc, aperiodic_auc, max_auc

type = 'max'; %periodic, aperiodic, max



% Plotting
figure
set(gcf,'renderer','painters');
scatter(good_outcome.(drs_to_compare),good_outcome.complete_auc,100,'filled')
hold on
scatter(bad_outcome.(drs_to_compare),bad_outcome.complete_auc,100,'filled')

correlat=corr(auc_tbl.(drs_to_compare),auc_tbl.complete_auc,'rows','complete','Type','Spearman');
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



