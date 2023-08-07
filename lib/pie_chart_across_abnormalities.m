%% iEEG%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scatter of abnormalities iEEG
clear all
close all

%%

cd /Users/c2056366/Documents/MATLAB/normativeBPToGMS/
parc=2;
%figdir='./figures/';
bands={'del';'the';'alp';'bet';'gam'};
load('lib/colourmaps_yw.mat');

%addpath(genpath('/Users/peter/Documents/MATLAB/altmany-export_fig-4c015d5'));
addpath(genpath('./scripts/'));
addpath(genpath('./lib/'));

resectedThresh=0.25;
%% preallocate space for MasterChannelTable
preallocate_table
%% UCLH
load_UCLH
%% RAM
load_RAM
%%
reorder_table_atlas
%%
 load_surface
%%
remove_invalid
%% Metadata
IDP={};ILAE1=[];
m=table(IDP,ILAE1);
clear metadata
load /Users/c2056366/Documents/MATLAB/normativeBPToGMS/data/metadata.mat
startlocs=size(m.IDP,1)+1;
endlocs=startlocs+(size(metadata,1)-1);
m.IDP(startlocs:endlocs)=convertStringsToChars(metadata.IDP); 
m.ILAE1(startlocs:endlocs)=metadata.ILAE1;
m.tle=metadata.istle;
m.etle=metadata.isfle;


metadata=m;
clear m startlocs endlocs
UCLsubjID=metadata.IDP;
nbUCLsubj=size(metadata,1);


%% Load fooof output
load_fooof_output_30hz

figdir='/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Figure3/';
load('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/UCLROIcontainsresected.mat')
%% load max abnormalities
load('/Users/c2056366/Documents/Brain_paper_data/Abnormalities/Max_z_score_aperiodic_fit.csv')
load('/Users/c2056366/Documents/Brain_paper_data/Abnormalities/Max_z_score_flat.csv')
load('/Users/c2056366/Documents/Brain_paper_data/Abnormalities/Max_z_score_fooof.csv')
load('/Users/c2056366/Documents/Brain_paper_data/Abnormalities/Max_z_score_slopes.csv')
load('/Users/c2056366/Documents/Brain_paper_data/Abnormalities/Max_z_score_original.csv')
load('/Users/c2056366/Documents/Brain_paper_data/Abnormalities/Abnormalities_4freqs_ieeg.mat')
%% load 30-47.5Hz slope abnormality

load('/Users/c2056366/Documents/Brain_paper_data/Abnormalities/Abnormalities_slope_ieeg_30_47_5.mat')
%% colormap and surface set up
cd('/Users/c2056366/Documents/iowaDATA/scripts/')
load('atlas.mat')
load('colourmaps_yw.mat')
T=atlas;
parc=2;
load_surface
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


%% Beeswarm plots
bad_patients=max_abn_patients(:,(metadata.ILAE1>2));
good_patients=abs(UCLROIdata_abnormality(:,(metadata.ILAE1<3)));
bad_patients_idx=selected_abn_patients(:,(metadata.ILAE1>2));
good_patients_idx=selected_abn_patients(:,(metadata.ILAE1<3));
bad_patient_res=UCLROIcontainsresected(:,(metadata.ILAE1>2));
good_patient_res=UCLROIcontainsresected(:,(metadata.ILAE1<3));
meta_good=metadata.IDP(metadata.ILAE1<3);
meta_bad=metadata.IDP(bad_idx);

FigH = figure('Position', get(0, 'Screensize'));
for i=1:34


    resected=find(good_patient_res(:,i)>resectedThresh);
    spared=find(good_patient_res(:,i)<=resectedThresh & good_patient_res(:,i)~=-1);
    if isempty(spared)
        continue
    end
    subplot(6,6,i)

    [~,~,~,drs]=perfcurve([zeros(size(resected));ones(size(spared))],good_patients([resected;spared],i),1);

    set(gcf,'position', [10 10 300 300])
    set(gcf,'renderer','painters');
    hold on
    [x,y]=UnivarScatter(padWithNans(good_patients(resected,i),good_patients(spared,i)),'Labels',{'resected','spared'},'MarkerFaceColor',cmaps.c);
    x(isnan(x))=[];
    y(isnan(y))=[];

    c=ceil((y./4)*256);
    c(isnan(c))=[];
    hold on
    scatter(x(:),y(:),ones(size(y(:),1),1)*200,cmaps.z_score_map(c(:),:),'filled')

    xlim([.5 2.5])
    ylim([0 4])
    caxis([0 4])

    title(meta_good{i},strcat(' D_R_S=',num2str(drs)))
    ylabel('max |z|')
    box off
end
f = gcf;
exportgraphics(f,'/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/iEEG_good_patients_abn_selected_beeswarm.png','Resolution',600)
close all
%% Collect the percentage of chosen abnomality per patient%colormap set up
colormap_idx=nan(5,3);
colormap_idx(1,:)=cmaps.delta_map(111,:);
colormap_idx(2,:)=cmaps.theta_map(111,:);
colormap_idx(3,:)=cmaps.alpha_map(111,:);
colormap_idx(4,:)=cmaps.beta_map(111,:);
map = brewermap(256,'Oranges'); 
map(1,:) = [1 1 1];
colormap_idx(5,:)=map(111,:);
dist_abn_patient=nan(63,5);

for i=1:63
    for z=1:5
        dist_abn_patient(i,z)= (sum(selected_abn_patients(:,i)==z)/sum(~isnan(selected_abn_patients(:,i))))*100;

    end
    subplot(8,8,i)
    pie2(dist_abn_patient(i,:))
    colormap(colormap_idx)
end




%% Outcome based
bad_idx=find((metadata.ILAE1>2));
good_idx=find((metadata.ILAE1<3));
meta_good=metadata.IDP(good_idx);
meta_bad=metadata.IDP(bad_idx);
bad_patients=selected_abn_patients(:,(metadata.ILAE1>2));
good_patients=selected_abn_patients(:,(metadata.ILAE1<3));
bad_patient_res=UCLROIcontainsresected(:,(metadata.ILAE1>2));
good_patient_res=UCLROIcontainsresected(:,(metadata.ILAE1<3));

%Implant based
tle_idx=find((metadata.tle==1));
etle_idx=find((metadata.etle==1));
meta_tle=metadata.IDP(tle_idx);
meta_etle=metadata.IDP(etle_idx);
tle_patients=selected_abn_patients(:,(metadata.tle==1));
etle_patients=selected_abn_patients(:,(metadata.etle==1));
tle_patient_res=UCLROIcontainsresected(:,(metadata.tle==1));
etle_patient_res=UCLROIcontainsresected(:,(metadata.etle==1));


FigH = figure('Position', get(0, 'Screensize'));
for i=1:29
   


resected=find(bad_patient_res(:,i)>resectedThresh);
spared=find(bad_patient_res(:,i)<=resectedThresh & bad_patient_res(:,i)~=-1);
if isempty(spared)
    continue
end
subplot(5,6,i)
subtitle(meta_bad{i})
set(gcf,'Color','w')
%set(gcf,'renderer','painters');
hold on
scatter3(T.xyz{parc}(spared,1),T.xyz{parc}(spared,2),T.xyz{parc}(spared,3),ones(size(spared,1),1)*50,[1 1 1],'filled');
scatter3(T.xyz{parc}(resected,1),T.xyz{parc}(resected,2),T.xyz{parc}(resected,3),ones(size(resected,1),1)*50,[0 0 0],'filled');
scatter3(T.xyz{parc}(:,1),T.xyz{parc}(:,2),T.xyz{parc}(:,3),ones(nbroi,1)*30,bad_patients(:,i),'filled');

colormap(colormap_idx)
cb=colorbar;
caxis([1 5])
set(cb,'FontSize',30)
set(cb,'YTick',0:1:4)
makeFigNice
hold on
resectedroilocs=[];
for i=1:numel(resected)
    resectedroilocs=[resectedroilocs;find(labelfaceL==resected(i))];
end
tmp=zeros(size(labelfaceL,1),1);
tmp(resectedroilocs)=1;
resectedroilocs=tmp;

 trisurf(lhpialface(find(resectedroilocs),:),lhpialvert(:,1),lhpialvert(:,2),lhpialvert(:,3)-140,'facecolor',[0.8 0.8 0.8],'facealpha',.2,'edgealpha',0)
  trisurf(lhpialface(find(~resectedroilocs),:),lhpialvert(:,1),lhpialvert(:,2),lhpialvert(:,3)-140,'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)
  trisurf(rhpialface,rhpialvert(:,1),rhpialvert(:,2),rhpialvert(:,3)-140,'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)

view([0 90]);
xlim([min(T.xyz{parc}(:))-20 max(T.xyz{parc}(:))+20]);
ylim([min(T.xyz{parc}(:))-20 max(T.xyz{parc}(:))+20]);
zlim([min(T.xyz{parc}(:))-120 max(T.xyz{parc}(:))+20]);
axis square
grid off
axis off
colorbar off
% title('max |z| across bands')
%         export_fig(strcat(figdir,'fig3brain',metadata.IDP{1},'.png'),'-png','-transparent','-m6');    



end

f = gcf;
exportgraphics(f,'/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/iEEG_bad_patients_abn_selected_ROI.png','Resolution',600)
close all
%% Patient based
bad_patients=selected_abn_patients(:,(metadata.ILAE1>2));
good_patients=selected_abn_patients(:,(metadata.ILAE1<3));

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

print(gcf, '-dpdf', '/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/iEEG_good_patients_abn_selected.pdf');



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

print(gcf, '-dpdf', '/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/iEEG_bad_patients_abn_selected.pdf');



%% Max calculation across tissue type

abnormality_dist_good_resected=nan(128,63);
abnormality_dist_good_spared=nan(128,63);
abnormality_dist_bad_resected=nan(128,63);
abnormality_dist_bad_spared=nan(128,63);

for i=1:63
    for z=1:128

        if (metadata.ILAE1(i)>2) & (Max_z_score_flat(z,i) > abs(Max_z_score_slopes(z,i))) & UCLROIcontainsresected(z,i)>resectedThresh
            abnormality_dist_bad_resected(z,i)=1;
        elseif (metadata.ILAE1(i)>2) & (Max_z_score_flat(z,i) < abs(Max_z_score_slopes(z,i))) & UCLROIcontainsresected(z,i)>resectedThresh
            abnormality_dist_bad_resected(z,i)=0;

        elseif (metadata.ILAE1(i)>2) & (Max_z_score_flat(z,i) > abs(Max_z_score_slopes(z,i))) & (UCLROIcontainsresected(z,i)<=resectedThresh & UCLROIcontainsresected(z,i)~=-1)
            abnormality_dist_bad_spared(z,i)=1;
        elseif (metadata.ILAE1(i)>2) & (Max_z_score_flat(z,i) < abs(Max_z_score_slopes(z,i))) & (UCLROIcontainsresected(z,i)<=resectedThresh & UCLROIcontainsresected(z,i)~=-1)
            abnormality_dist_bad_spared(z,i)=0;

        elseif (metadata.ILAE1(i)<3) & (Max_z_score_flat(z,i) > abs(Max_z_score_slopes(z,i))) & UCLROIcontainsresected(z,i)>resectedThresh
            abnormality_dist_good_resected(z,i)=1;
        elseif (metadata.ILAE1(i)<3) & (Max_z_score_flat(z,i) < abs(Max_z_score_slopes(z,i))) & UCLROIcontainsresected(z,i)>resectedThresh
            abnormality_dist_good_resected(z,i)=0;
            
        elseif (metadata.ILAE1(i)<3) & (Max_z_score_flat(z,i) > abs(Max_z_score_slopes(z,i))) & (UCLROIcontainsresected(z,i)<=resectedThresh & UCLROIcontainsresected(z,i)~=-1)
            abnormality_dist_good_spared(z,i)=1;
        elseif (metadata.ILAE1(i)<3) & (Max_z_score_flat(z,i) < abs(Max_z_score_slopes(z,i))) & (UCLROIcontainsresected(z,i)<=resectedThresh & UCLROIcontainsresected(z,i)~=-1)
            abnormality_dist_good_spared(z,i)=0;
        end
    end
end



%% concatonate them and labels gen and pie plot
labels={'Periodic Bandpower','Aperiodic exponent'};
X=nan(1,2);

%%
collapsed_good_resected=cat(1,abnormality_dist_good_resected(:));
collapsed_good_resected=rmmissing(collapsed_good_resected);

X(1)=sum(collapsed_good_resected)/size(collapsed_good_resected,1);
X(2)=(size(collapsed_good_resected,1)-sum(collapsed_good_resected))/size(collapsed_good_resected,1);
X(1)=nanmean(rati_tabs(:,1));
X(2)=1-nanmean(rati_tabs(:,1));
subplot(2,2,1)


pie(X)
title('Max z chosen in ILAE 1,2 resected cases')

%%
collapsed_bad_resected=cat(1,abnormality_dist_bad_resected(:));
collapsed_bad_resected=rmmissing(collapsed_bad_resected);

X(1)=sum(collapsed_bad_resected)/size(collapsed_bad_resected,1);
X(2)=(size(collapsed_bad_resected,1)-sum(collapsed_bad_resected))/size(collapsed_bad_resected,1);
X(1)=nanmean(rati_tabs(:,2));
X(2)=1-nanmean(rati_tabs(:,2));
subplot(2,2,2)

pie(X)
title('Max z chosen in ILAE 3+ resected cases')

%%
collapsed_good_spared=cat(1,abnormality_dist_good_spared(:));
collapsed_good_spared=rmmissing(collapsed_good_spared);

X(1)=sum(collapsed_good_spared)/size(collapsed_good_spared,1);
X(2)=(size(collapsed_good_spared,1)-sum(collapsed_good_spared))/size(collapsed_good_spared,1);
X(1)=nanmean(rati_tabs(:,3));
X(2)=1-nanmean(rati_tabs(:,3));
subplot(2,2,3)

pie(X)
title('Max z chosen in ILAE 1,2 spared cases')

%%
collapsed_bad_spared=cat(1,abnormality_dist_bad_spared(:));
collapsed_bad_spared=rmmissing(collapsed_bad_spared);

X(1)=sum(collapsed_bad_spared)/size(collapsed_bad_spared,1);
X(2)=(size(collapsed_bad_spared,1)-sum(collapsed_bad_spared))/size(collapsed_bad_spared,1);
X(1)=nanmean(rati_tabs(:,4));
X(2)=1-nanmean(rati_tabs(:,4));
subplot(2,2,4)

pie(X)
title('Max z chosen in ILAE 3+ spared cases')


% Create legend
lgd = legend(labels);
lgd.Layout.Tile = 'northeastoutside';

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);
set(h,'renderer','painters');
print(gcf, '-dpdf', '/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Figure3/iEEG_Piechart_dist_of_max_chosen.pdf');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MEG%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
figdir='/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Figure4';

%% 

MEG_auc=readtable('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/MEG_drs/MEG_auc_cs_mixed.csv');
metadata_full=readtable('/Users/c2056366/Documents/MEG_Tom/meta_meg.csv');
load('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Max_abnormality_across_comps/MEG_abnormality_matrices/MEG_resection.mat');
neo_resected=neo_resected';
MEG_auc = sortrows(MEG_auc,"ID","ascend");
excluded_pats=metadata_full.ID(~ismember(metadata_full.ID,MEG_auc.ID));
metadata_full(ismember(metadata_full.ID,excluded_pats),:)=[];

ETLE_idx=find(strcmp(metadata_full.isNeocortical,'Yes'));
MEG_auc=MEG_auc(ETLE_idx,:);

metadata=table();
metadata.ID=MEG_auc.ID;
metadata.ILAE1=MEG_auc.outcome;
metadata.ILAE1(strcmp(metadata.ILAE1,'ILAE1'))={1};
metadata.ILAE1(strcmp(metadata.ILAE1,'ILAE2+'))={3};
metadata.ILAE1=cell2mat(metadata.ILAE1);

%% Load abnromality matrices
load('/Users/c2056366/Documents/MEG_Tom/ucl_bandpower_epoch1_csaba.mat')
IDs=bandpower_results(1).subjID;

ELTE_IDs=IDs(ismember(IDs,metadata.ID));
data_etle_idx=find(ismember(IDs,ELTE_IDs));

load('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Max_abnormality_across_comps/MEG_abnormality_matrices/MEG_max_abnormalities_flat.csv')
load('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Max_abnormality_across_comps/MEG_abnormality_matrices/MEG_max_abnormalities_slope.csv')
load('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Max_abnormality_across_comps/MEG_abnormality_matrices/MEG_max_abnormalities.csv')
load('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Max_abnormality_across_comps/MEG_abnormality_matrices/MEG_abnormalities_freqs.mat');
MEG_freq_abnormalities=foo;
MEG_max_abnormalities=MEG_max_abnormalities(:,data_etle_idx);
MEG_max_abnormalities_flat=MEG_max_abnormalities_flat(:,data_etle_idx);
MEG_max_abnormalities_slope=MEG_max_abnormalities_slope(:,data_etle_idx);

MEG_abnormality=table();
MEG_abnormality.ID=ELTE_IDs;
MEG_abnormality.Max_abnormality=MEG_max_abnormalities';
MEG_abnormality.Slope_abnormality=MEG_max_abnormalities_slope';
MEG_abnormality.Flat_abnormality=MEG_max_abnormalities_flat';
MEG_abnormality=sortrows(MEG_abnormality,"ID","ascend");
MEG_abnormality.ILAE=metadata.ILAE1;

%% Max calculation across each freqs and aperiodic
selected_abn_patients=nan(114,33);

for z=1:33
    for i=1:114
        if isnan(MEG_freq_abnormalities(i,1,z))
            continue
        end
           


     
            mx_indi=[abs(MEG_freq_abnormalities(i,1,z)), abs(MEG_freq_abnormalities(i,2,z)),abs(MEG_freq_abnormalities(i,3,z)),abs(MEG_freq_abnormalities(i,4,z)),abs(MEG_max_abnormalities_slope(i,z))];
            [M,I] =max(mx_indi);
           selected_abn_patients(i,z)=I;
        
    end
    
end
%% %% Beeswarm plots
load('/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Max_abnormality_across_comps/MEG_abnormality_matrices/MEG_max_abnormalities_slope_30_47_5.csv')
load('/Users/c2056366/Documents/MEG_Tom/ucl_bandpower_epoch1_csaba.mat')
IDs=bandpower_results(1).subjID;
resectedThresh=0.25;
ELTE_IDs=IDs(ismember(IDs,metadata.ID));
data_etle_idx=find(ismember(IDs,ELTE_IDs));

MEG_max_abnormalities=MEG_max_abnormalities_slope_30_47_5(:,data_etle_idx);

MEG_abnormality=table();
MEG_abnormality.ID=ELTE_IDs;
MEG_abnormality.Max_abnormality=MEG_max_abnormalities';
MEG_abnormality=sortrows(MEG_abnormality,"ID","ascend");
MEG_abnormality.ILAE=metadata.ILAE1;

good_patients=abs(MEG_abnormality.Max_abnormality((MEG_abnormality.ILAE<3),:));
good_patient_res=neo_resected((MEG_abnormality.ILAE<3),:);
meta_good=MEG_abnormality.ID(MEG_abnormality.ILAE<3);
%%
FigH = figure('Position', get(0, 'Screensize'));
for i=1:12


    resected=find(good_patient_res(i,:)>resectedThresh);
    spared=find(good_patient_res(i,:)<=resectedThresh);
    if isempty(spared) || isempty(resected) 
        continue
    end
    subplot(3,4,i)

    [~,~,~,drs]=perfcurve([zeros(size(resected)),ones(size(spared))],good_patients(i,[resected,spared]),1);

    set(gcf,'position', [10 10 300 300])
    set(gcf,'renderer','painters');
    hold on
    [x,y]=UnivarScatter(padWithNans(good_patients(i,resected),good_patients(i,spared)),'Labels',{'resected','spared'},'MarkerFaceColor',cmaps.c);
    x(isnan(x))=[];
    y(isnan(y))=[];

    c=ceil((y./7)*256);
    c(isnan(c))=[];
    hold on
    scatter(x(:),y(:),ones(size(y(:),1),1)*200,cmaps.z_score_map(c(:),:),'filled')

    xlim([.5 2.5])
    ylim([0 4])
    clim([0 4])

    title(num2str(meta_good(i)),strcat(' D_R_S=',num2str(drs)))
    ylabel('max |z|')
    box off
end
%%
f = gcf;
exportgraphics(f,'/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/iEEG_good_patients_abn_selected_beeswarm.png','Resolution',600)
close all

%% ROI based
bad_idx=find((metadata.ILAE1>2));
good_idx=find((metadata.ILAE1<3));
bad_patients=selected_abn_patients(:,(metadata.ILAE1>2));
good_patients=selected_abn_patients(:,(metadata.ILAE1<3));

bad_patient_res=neo_resected((metadata.ILAE1>2),:);
good_patient_res=neo_resected((metadata.ILAE1<3),:);
meta_good=metadata.ID(good_idx);
meta_bad=metadata.ID(bad_idx);

resectedThresh=.25;
cd('/Users/c2056366/Documents/iowaDATA/scripts/')
load('atlas.mat')
load('colourmaps_yw.mat')
T=atlas;
parc=2;
load_surface

T.xyz{2,1}=T.xyz{2,1}(15:end,:);

%colormap set up
colormap_idx=nan(5,3);
colormap_idx(1,:)=cmaps.delta_map(111,:);
colormap_idx(2,:)=cmaps.theta_map(111,:);
colormap_idx(3,:)=cmaps.alpha_map(111,:);
colormap_idx(4,:)=cmaps.beta_map(111,:);
map = brewermap(256,'Oranges'); 
map(1,:) = [1 1 1];
colormap_idx(5,:)=map(111,:);



FigH = figure('Position', get(0, 'Screensize'));
for i=1:12
   


resected=find(good_patient_res(:,i)>resectedThresh);
spared=find(good_patient_res(:,i)<=resectedThresh) ;


subplot(3,7,i)
title(num2str(meta_good(i)))
set(gcf,'Color','w')
%set(gcf,'renderer','painters');
hold on
scatter3(T.xyz{parc}(spared,1),T.xyz{parc}(spared,2),T.xyz{parc}(spared,3),ones(size(spared,1),1)*70,[1 1 1],'filled');
scatter3(T.xyz{parc}(resected,1),T.xyz{parc}(resected,2),T.xyz{parc}(resected,3),ones(size(resected,1),1)*70,[0 0 0],'filled');
scatter3(T.xyz{parc}(:,1),T.xyz{parc}(:,2),T.xyz{parc}(:,3),ones(114,1)*40,good_patients(:,i),'filled');

colormap(colormap_idx)
cb=colorbar;
caxis([1 5])
set(cb,'FontSize',30)
set(cb,'YTick',0:1:4)
makeFigNice
hold on
resectedroilocs=[];
for i=1:numel(resected)
    resectedroilocs=[resectedroilocs;find(labelfaceL==resected(i))];
end
tmp=zeros(size(labelfaceL,1),1);
tmp(resectedroilocs)=1;
resectedroilocs=tmp;

 trisurf(lhpialface(find(resectedroilocs),:),lhpialvert(:,1),lhpialvert(:,2),lhpialvert(:,3)-140,'facecolor',[0.8 0.8 0.8],'facealpha',.2,'edgealpha',0)
  trisurf(lhpialface(find(~resectedroilocs),:),lhpialvert(:,1),lhpialvert(:,2),lhpialvert(:,3)-140,'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)
  trisurf(rhpialface,rhpialvert(:,1),rhpialvert(:,2),rhpialvert(:,3)-140,'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)

view([0 90]);
xlim([min(T.xyz{parc}(:))-20 max(T.xyz{parc}(:))+20]);
ylim([min(T.xyz{parc}(:))-20 max(T.xyz{parc}(:))+20]);
zlim([min(T.xyz{parc}(:))-120 max(T.xyz{parc}(:))+20]);
axis square
grid off
axis off
colorbar off
 title('max |z| across bands')
%         export_fig(strcat(figdir,'fig3brain',metadata.IDP{1},'.png'),'-png','-transparent','-m6');    

title(num2str(meta_good(i)))

end

f = gcf;
exportgraphics(f,'/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/MEG_good_patients_abn_selected_ROI.png','Resolution',600)
close all
%% Patient based
bad_idx=find((metadata.ILAE1>2));
good_idx=find((metadata.ILAE1<3));
bad_patients=selected_abn_patients(:,(metadata.ILAE1>2));
good_patients=selected_abn_patients(:,(metadata.ILAE1<3));


load('colourmaps_yw.mat')

%colormap set up
colormap_idx=nan(5,3);
colormap_idx(1,:)=cmaps.delta_map(111,:);
colormap_idx(2,:)=cmaps.theta_map(111,:);
colormap_idx(3,:)=cmaps.alpha_map(111,:);
colormap_idx(4,:)=cmaps.beta_map(111,:);
map = brewermap(256,'Oranges'); 
map(1,:) = [1 1 1];
colormap_idx(5,:)=map(111,:);
meta_good=metadata.ID(good_idx);
meta_bad=metadata.ID(bad_idx);


%%
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

print(gcf, '-dpdf', '/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/MEG_good_patients_abn_selected.pdf');



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

print(gcf, '-dpdf', '/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Supplementary/MEG_bad_patients_abn_selected.pdf');



%% Max calculation across tissue type

resectedThresh=0.25;
abnormality_dist_good_resected=nan(33,114);
abnormality_dist_good_spared=nan(33,114);
abnormality_dist_bad_resected=nan(33,114);
abnormality_dist_bad_spared=nan(33,114);

for z=1:33
    for i=1:114

        if (MEG_abnormality.ILAE(z)>2) & (MEG_abnormality.Flat_abnormality(z,i) > abs(MEG_abnormality.Slope_abnormality(z,i))) & neo_resected(z,i)>resectedThresh
            abnormality_dist_bad_resected(z,i)=1;
        elseif (MEG_abnormality.ILAE(z)>2) & (MEG_abnormality.Flat_abnormality(z,i) < abs(MEG_abnormality.Slope_abnormality(z,i))) & neo_resected(z,i)>resectedThresh
            abnormality_dist_bad_resected(z,i)=0;

        elseif (MEG_abnormality.ILAE(z)>2) & (MEG_abnormality.Flat_abnormality(z,i) > abs(MEG_abnormality.Slope_abnormality(z,i))) & (neo_resected(z,i)<=resectedThresh )
            abnormality_dist_bad_spared(z,i)=1;
        elseif (MEG_abnormality.ILAE(z)>2) & (MEG_abnormality.Flat_abnormality(z,i) < abs(MEG_abnormality.Slope_abnormality(z,i))) & (neo_resected(z,i)<=resectedThresh )
            abnormality_dist_bad_spared(z,i)=0;

        elseif (MEG_abnormality.ILAE(z)<3) & (MEG_abnormality.Flat_abnormality(z,i) > abs(MEG_abnormality.Slope_abnormality(z,i))) & neo_resected(z,i)>resectedThresh
            abnormality_dist_good_resected(z,i)=1;
        elseif (MEG_abnormality.ILAE(z)<3) & (MEG_abnormality.Flat_abnormality(z,i) < abs(MEG_abnormality.Slope_abnormality(z,i))) & neo_resected(z,i)>resectedThresh
            abnormality_dist_good_resected(z,i)=0;
            
        elseif (MEG_abnormality.ILAE(z)<3) & (MEG_abnormality.Flat_abnormality(z,i) > abs(MEG_abnormality.Slope_abnormality(z,i))) & (neo_resected(z,i)<=resectedThresh )
            abnormality_dist_good_spared(z,i)=1;
        elseif (MEG_abnormality.ILAE(z)<3) & (MEG_abnormality.Flat_abnormality(z,i) < abs(MEG_abnormality.Slope_abnormality(z,i))) & (neo_resected(z,i)<=resectedThresh )
            abnormality_dist_good_spared(z,i)=0;
        end
    end
end

%% Ratio of chosen per patient
bad_resected=cell(33,1);
bad_spared=cell(33,1);
good_resected=cell(33,1);
good_spared=cell(33,1);


bad_resected_idx=(~isnan(abnormality_dist_bad_resected));
bad_spared_idx=(~isnan(abnormality_dist_bad_spared));
good_resected_idx=(~isnan(abnormality_dist_good_resected));
good_spared_idx=(~isnan(abnormality_dist_good_spared));

ratio_bad_resected=nan(33,1);
ratio_bad_spared=nan(33,1);
ratio_good_resected=nan(33,1);
ratio_good_spared=nan(33,1);

for i=1:33
    if sum(bad_resected_idx(:,i))>0
        bad_resected{i}=abnormality_dist_bad_resected(bad_resected_idx(:,i)>0,i);
    end
    if sum(bad_spared_idx(:,i))>0
        bad_spared{i}=abnormality_dist_bad_spared(bad_spared_idx(:,i)>0,i);
    end
    if sum(good_resected_idx(:,i))>0
        good_resected{i}=abnormality_dist_good_resected(good_resected_idx(:,i)>0,i);
    end
    if sum(good_spared_idx(:,i))>0
        good_spared{i}=abnormality_dist_good_spared(good_spared_idx(:,i)>0,i);
    end

    if ~isempty(bad_resected{i})
        ratio_bad_resected(i)=sum(bad_resected{i})/size(bad_resected{i},1);
    end
    if ~isempty(bad_spared{i})
        ratio_bad_spared(i)=sum(bad_spared{i})/size(bad_spared{i},1);
    end
    if ~isempty(good_resected{i})
        ratio_good_resected(i)=sum(good_resected{i})/size(good_resected{i},1);
    end
    if ~isempty(good_spared{i})
        ratio_good_spared(i)=sum(good_spared{i})/size(good_spared{i},1);
    end
end

%% corr
rati_tabs = zeros(size(ratio_good_spared,1),4);
rati_tabs(:,1) =ratio_bad_resected;
rati_tabs(:,2) = ratio_bad_spared;
rati_tabs(:,3) = ratio_good_resected;
rati_tabs(:,4) = ratio_good_spared;

% Replace upper triangle with NaNs
c=redblueTecplot;
h = heatmap(rati_tabs,'MissingDataColor','w','Colormap',c);
labels = ["Bad outcome Resected ROIs","Bad outcome Spared ROIs","Good outcome Resected ROIs","Good outcome Spared ROIs"];
h.XDisplayLabels = labels;
title('Proportion of iEEG ROIs where Periodic abnormality was chosen')
ylabel('Patients')
h=gcf;
%set(h,'PaperOrientation','landscape');
%set(h,'PaperPosition', [1 1 28 19]);
set(h,'renderer','painters');


%% concatonate them and labels gen and pie plot
labels={'Periodic Bandpower','Aperiodic exponent'};
X=nan(1,2);

%%
collapsed_good_resected=cat(1,abnormality_dist_good_resected(:));
collapsed_good_resected=rmmissing(collapsed_good_resected);

X(1)=sum(collapsed_good_resected)/size(collapsed_good_resected,1);
X(2)=(size(collapsed_good_resected,1)-sum(collapsed_good_resected))/size(collapsed_good_resected,1);
X(1)=nanmean(rati_tabs(:,1));
X(2)=1-nanmean(rati_tabs(:,1));
subplot(2,2,1)


pie(X)
title('Max z chosen in ILAE 1,2 resected cases')

%%
collapsed_bad_resected=cat(1,abnormality_dist_bad_resected(:));
collapsed_bad_resected=rmmissing(collapsed_bad_resected);

X(1)=sum(collapsed_bad_resected)/size(collapsed_bad_resected,1);
X(2)=(size(collapsed_bad_resected,1)-sum(collapsed_bad_resected))/size(collapsed_bad_resected,1);
X(1)=nanmean(rati_tabs(:,2));
X(2)=1-nanmean(rati_tabs(:,2));
subplot(2,2,2)

pie(X)
title('Max z chosen in ILAE 3+ resected cases')

%%
collapsed_good_spared=cat(1,abnormality_dist_good_spared(:));
collapsed_good_spared=rmmissing(collapsed_good_spared);

X(1)=sum(collapsed_good_spared)/size(collapsed_good_spared,1);
X(2)=(size(collapsed_good_spared,1)-sum(collapsed_good_spared))/size(collapsed_good_spared,1);
X(1)=nanmean(rati_tabs(:,3));
X(2)=1-nanmean(rati_tabs(:,3));
subplot(2,2,3)

pie(X)
title('Max z chosen in ILAE 1,2 spared cases')

%%
collapsed_bad_spared=cat(1,abnormality_dist_bad_spared(:));
collapsed_bad_spared=rmmissing(collapsed_bad_spared);

X(1)=sum(collapsed_bad_spared)/size(collapsed_bad_spared,1);
X(2)=(size(collapsed_bad_spared,1)-sum(collapsed_bad_spared))/size(collapsed_bad_spared,1);
X(1)=nanmean(rati_tabs(:,4));
X(2)=1-nanmean(rati_tabs(:,4));
subplot(2,2,4)

pie(X)
title('Max z chosen in ILAE 3+ spared cases')


% Create legend
lgd = legend(labels);
lgd.Layout.Tile = 'northeastoutside';

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);
set(h,'renderer','painters');
print(gcf, '-dpdf', '/Users/c2056366/Documents/iowaDATA/fooof_output/For paper/Figures/Figure4/MEG_Piechart_dist_of_max_chosen.pdf');