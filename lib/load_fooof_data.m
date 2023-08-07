function [normative_table,patient_table]=load_fooof_data(datapath)

% Load all fooof outputs for abnormality calculation
%% Load aperiodic components
aperiodic_RAM=readtable([datapath,'RAM_aperiodic_30hz.csv']);
aperiodic_UCLH=readtable([datapath,'UCLH_aperiodic_30hz.csv']);
%% Load flattened PSD
pxx_RAM_flattened_psd_30hz=load([datapath,'pxx_RAM_flattened_psd_30hz.csv']);
pxx_UCLH_flattened_psd_30hz=load([datapath,'pxx_UCLH_flattened_psd_30hz.csv']);

%% unifiy necessary columns
%RAM
normative_table=aperiodic_RAM;
normative_table.flattened_psd=pxx_RAM_flattened_psd_30hz;
normative_table.DataSource(:,:)="RAM";
%UCLH
patient_table=aperiodic_UCLH;
patient_table.flattened_psd=pxx_UCLH_flattened_psd_30hz;
patient_table.DataSource(:,:)="UCLH";
end

