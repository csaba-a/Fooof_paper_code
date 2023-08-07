%%%% Figure 2 %%%%%
clear all
close all

%% paths and features
define_paths

%% Rerun this section when assign new example patient
% Example patient
example_patient="1216";%1216 --good outcome ; 910 -- bad outcome


% Load and set up metadata and select example patient
[metadata,UCLsubjID,~]=load_metadata(metadata_path);

%select out the right example patient to analyse and plot
metadata(find(metadata.IDP~=example_patient),:)=[];
UCLsubjID(find(UCLsubjID~=example_patient),:)=[];

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

type_power_spectrum='max'; %complete, periodic, aperiodic, max


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
        abnormalities_tbl.complete_abn=UCLROIdata_abnormality_complete;

        % Plot abnormality brain plots
        [resected,spared]=plot_abnormalities(UCLROIdata_abnormality_complete,UCLROIcontainsresected,resectedThresh,parc);
        % Plot DRS plots
        drs_plot_indiv(resected,spared,UCLROIdata_abnormality_complete);

    case 'periodic' %Periodic
        [UCLROIdata_abnormality_periodic,UCLROIcontainsresected]=calc_abnormality(normative_table, patient_table,UCLsubjID,n_bands,resectedThresh, parc);
        abnormalities_tbl.periodic_abn=UCLROIdata_abnormality_periodic;

        % Plot abnormality brain plots
        [resected,spared]=plot_abnormalities(UCLROIdata_abnormality_periodic,UCLROIcontainsresected,resectedThresh,parc);
        % Plot DRS plots
        drs_plot_indiv(resected,spared,UCLROIdata_abnormality_periodic);

        
    case 'aperiodic' %Aperiodic
        [UCLROIdata_abnormality_aperiodic,UCLROIcontainsresected]=calc_abnormality_aperiodic(normative_table, patient_table,UCLsubjID,resectedThresh, parc);
        abnormalities_tbl.aperiodic_abn=UCLROIdata_abnormality_aperiodic;

        % Plot abnormality brain plots
        [resected,spared]=plot_abnormalities(UCLROIdata_abnormality_aperiodic,UCLROIcontainsresected,resectedThresh,parc);
        % Plot DRS plots
        drs_plot_indiv(resected,spared,UCLROIdata_abnormality_aperiodic);

       

    case 'max' %Max abnormality across periodic and aperiodic
        for z=1:length(UCLsubjID)
            UCLROIdata_abnormality_max(:,z)=max(abs(UCLROIdata_abnormality_periodic(:,z)),abs(UCLROIdata_abnormality_aperiodic(:,z)));
        end
            abnormalities_tbl.max_abn=UCLROIdata_abnormality_max;

            % Plot abnormality brain plots
            [resected,spared]=plot_abnormalities(UCLROIdata_abnormality_max,UCLROIcontainsresected,resectedThresh,parc);
            % Plot DRS plots
            drs_plot_indiv(resected,spared,UCLROIdata_abnormality_max);

           
        
end

% Save plots
save_figs(example_patient,type_power_spectrum,analysis_location,figdir_2)
%% Compare abnormality distributions with correlation plots %%
%Run after generated all the abnormality matrices

% Define the type of abnormality you want to compare to Complete Bp
% abnormality
abnormality_to_compare = 'aperiodic_abn'; %abnormalities_tbl.aperiodic_abn or abnormalities_tbl.max_abn
type = 'aperiodic'; %periodic, aperiodic, max
% Plotting

figure
set(gcf,'renderer','painters');
scatter(abnormalities_tbl.(abnormality_to_compare)(spared),abnormalities_tbl.complete_abn(spared),100,'filled')
hold on
scatter(abnormalities_tbl.(abnormality_to_compare)(resected),abnormalities_tbl.complete_abn(resected),100,'filled')

correlat=corr(abnormalities_tbl.(abnormality_to_compare),abnormalities_tbl.complete_abn,'rows','complete','Type','Spearman');
legend({'Spared','Resected'},'Location','northwest')
xlabel(['Max z based on ', type, ' PSD'])
ylabel('Max z score based on complete PSD')
set(gca, 'FontSize', 16)
title(['r=',num2str(round(correlat,2))])
%axis equal
set(gca,'XTick',[])
set(gca,'YTick',[])

%Save figure
saveas(gca, fullfile(analysis_location,figdir_2,strcat(example_patient,'_scatter_Plot_',type,'_vs_complete_abnormality', '.pdf'))); % specify the full path
close all
