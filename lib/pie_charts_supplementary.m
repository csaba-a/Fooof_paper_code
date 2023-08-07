clear all
close all
%% paths and features
addpath(genpath('/Users/c2056366/Documents/Fooof_paper_material')) %custom scripts
original_datapath='/Users/c2056366/Documents/Fooof_paper_material/iEEG_data/'; %iEEG data
fooof_output_path='/Users/c2056366/Documents/Fooof_paper_material/fooof_output/';%fooof output
metadata_path='/Users/c2056366/Documents/Fooof_paper_material/iEEG_data/data/metadata.mat'; %corresponding metadata
figdir='/Users/c2056366/Documents/Fooof_paper_material/Figures/Supplementary/';% figure save
iEEG_abn_path='/Users/c2056366/Documents/Fooof_paper_material/iEEG Abnormalities/'; % iEEG abnormalities path
%%
addpath(genpath('/Users/c2056366/Documents/fieldtrip-20221005')) %Fieldtrip
%%
% define features
fi=1:0.5:30;
parc=2;

% load colorscheme
load('lib/colourmaps_yw.mat');
% Resection threshold
resectedThresh=0.25;

% Load fooof output and set up path and periodic/original data
[normative_table,patient_table]=load_fooof_data(fooof_output_path);

% Load original data

[ATLAS, labelfaceL,labelfaceR,lhpialface,lhpialvert,lobes,MasterChannelTable,nbroi,rhpialface,rhpialvert,T]=prep_original_iEEG_data(original_datapath,parc);
atlas=T;

% Load and set up metadata and select example patient
[metadata,UCLsubjID,~]=load_metadata(metadata_path);
%% Load resection matrix

load([iEEG_abn_path,'resections.mat'])
%% Load abnormalities
load([iEEG_abn_path,'Max_z_score_aperiodic_fit.csv'])
load([iEEG_abn_path,'Max_z_score_flat.csv'])
load([iEEG_abn_path,'Max_z_score_fooof.csv'])
load([iEEG_abn_path,'Max_z_score_slopes.csv'])
load([iEEG_abn_path,'Max_z_score_original.csv'])
load([iEEG_abn_path,'Abnormalities_4freqs_ieeg.mat'])

%% colormap set up
colormap_idx=nan(5,3);
colormap_idx(1,:)=cmaps.delta_map(111,:);
colormap_idx(2,:)=cmaps.theta_map(111,:);
colormap_idx(3,:)=cmaps.alpha_map(111,:);
colormap_idx(4,:)=cmaps.beta_map(111,:);
map = brewermap(256,'Oranges');
map(1,:) = [1 1 1];
colormap_idx(5,:)=map(111,:);

%% Get the max abnormality
max_abn_patients=nan(128,63);
selected_abn_patients=nan(128,63);
for z=1:63
    for i=1:128
        if isnan(UCLROIdata_abnormality(i,1,z))
            continue
        end
        mx_indi=[abs(UCLROIdata_abnormality(i,1,z)), abs(UCLROIdata_abnormality(i,2,z)),abs(UCLROIdata_abnormality(i,3,z)),abs(UCLROIdata_abnormality(i,4,z)),abs(Max_z_score_aperiodic_fit(i,z))];
        [M,I] =max(mx_indi);
        max_abn_patients(i,z)=M;
        selected_abn_patients(i,z)=I;

    end

end

%% Calcualte and plot pie chart of chosen frequency bands or aperiodic exponent
% Outcome based iEEG
%% Patient based
bad_patients=selected_abn_patients(:,(metadata.ILAE1>2));
good_patients=selected_abn_patients(:,(metadata.ILAE1<3));
bad_idx=find((metadata.ILAE1>2));
good_idx=find((metadata.ILAE1<3));
meta_good=metadata.IDP(good_idx);
meta_bad=metadata.IDP(bad_idx);


for i=1:size(good_patients,2)
    tbl=tabulate(good_patients(:,i));
    subplot(5,7,i)
    p=pie(tbl(:,3));
    colormap(colormap_idx)
    delete(findobj(p,'Type','text'))%option 2
    subtitle(meta_good{i})
end

h=gcf;
set(h,'renderer','painters');
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 29 19]);

print(gcf, '-dpdf', [figdir,'iEEG_good_patients_abn_selected.pdf']);



for i=1:size(bad_patients,2)
    tbl=tabulate(bad_patients(:,i));
    subplot(5,6,i)
    p=pie(tbl(:,3));
    colormap(colormap_idx)
    delete(findobj(p,'Type','text'))%option 2
    subtitle(meta_bad{i})
end

h=gcf;
set(h,'renderer','painters');
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 29 19]);

print(gcf, '-dpdf', [figdir,'iEEG_bad_patients_abn_selected.pdf']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MEG%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except figdir colormap_idx
close all
%% Load paths
addpath(genpath('/Users/c2056366/Documents/Fooof_paper_material')) %custom scripts
addpath(genpath('/Users/c2056366/Documents/fieldtrip-20221005')) %Fieldtrip
figdir='/Users/c2056366/Documents/Fooof_paper_material/Figures/Figure4/';% figure save
meg_data_path = '/Users/c2056366/Documents/Fooof_paper_material/MEG_data/';
MEG_abn_path='/Users/c2056366/Documents/Fooof_paper_material/MEG Abnormalities/'; % iEEG abnormalities path

%% Load Drs values and metadata

%Drs values
[MEG_auc]=load_meg_aucs(meg_data_path);

%Metadata
[metadata,MEG_auc]=loadmetadata_MEG(meg_data_path,MEG_auc);
%% Load IDs, abnormality matrices
load([MEG_abn_path,'MEG_patient_IDs.mat']);

ELTE_IDs=IDs(ismember(IDs,metadata.ID));
data_etle_idx=find(ismember(IDs,ELTE_IDs));

load([MEG_abn_path,'MEG_max_abnormalities_flat.csv'])
load([MEG_abn_path,'MEG_max_abnormalities_slope.csv'])
load([MEG_abn_path,'MEG_max_abnormalities.csv'])
load([MEG_abn_path,'MEG_abnormalities_freqs.mat']);

MEG_freq_abnormalities=foo;
MEG_max_abnormalities=MEG_max_abnormalities(:,data_etle_idx);
MEG_max_abnormalities_flat=MEG_max_abnormalities_flat(:,data_etle_idx);
MEG_max_abnormalities_slope=MEG_max_abnormalities_slope(:,data_etle_idx);

%% Assamble abnormality table for plotting
MEG_abnormality=table();
MEG_abnormality.ID=ELTE_IDs;
MEG_abnormality.Max_abnormality=MEG_max_abnormalities';
MEG_abnormality.Slope_abnormality=MEG_max_abnormalities_slope';
MEG_abnormality.Flat_abnormality=MEG_max_abnormalities_flat';
MEG_abnormality=sortrows(MEG_abnormality,"ID","ascend");
MEG_abnormality.ILAE=metadata.ILAE1;

%% Max calculation across each freqs and aperiodic
selected_abn_patients=nan(114,size(MEG_abnormality,1));

for z=1:size(MEG_abnormality,1)
    for i=1:114
        if isnan(MEG_freq_abnormalities(i,1,z))
            continue
        end
        mx_indi=[abs(MEG_freq_abnormalities(i,1,z)), abs(MEG_freq_abnormalities(i,2,z)),abs(MEG_freq_abnormalities(i,3,z)),abs(MEG_freq_abnormalities(i,4,z)),abs(MEG_max_abnormalities_slope(i,z))];
        [M,I] =max(mx_indi);
        selected_abn_patients(i,z)=I;

    end

end
%% Separate patient data and metadata based on outcome
bad_idx=find((metadata.ILAE1>2));
good_idx=find((metadata.ILAE1<3));
bad_patients=selected_abn_patients(:,(metadata.ILAE1>2));
good_patients=selected_abn_patients(:,(metadata.ILAE1<3));
meta_good=metadata.ID(good_idx);
meta_bad=metadata.ID(bad_idx);


%% Plotting pie charts for MEG

% Good Patients
for i=1:size(good_patients,2)
    tbl=tabulate(good_patients(:,i));
    subplot(2,6,i)
    p=pie(tbl(:,3));
    colormap(colormap_idx)
    delete(findobj(p,'Type','text'))%option 2
    subtitle(meta_good(i))
end

h=gcf;
set(h,'renderer','painters');
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 29 19]);

print(gcf, '-dpdf', [figdir,'MEG_good_patients_abn_selected.pdf']);


% Bad patients
for i=1:size(bad_patients,2)
    tbl=tabulate(bad_patients(:,i));
    subplot(3,7,i)
    p=pie(tbl(:,3));
    colormap(colormap_idx)
    delete(findobj(p,'Type','text'))%option 2
    subtitle(meta_bad(i))
end

h=gcf;
set(h,'renderer','painters');
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 29 19]);

print(gcf, '-dpdf', [figdir,'MEG_bad_patients_abn_selected.pdf']);

