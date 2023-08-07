function [UCLROIdata_abnormality_final,UCLROIcontainsresected]=calc_abnormality(normative_table, patient_table,UCLsubjID,nbands,resectedThresh, parc)
%% CALC_ABNORMALITY
% Generates the normative map from control cohort via 
% taking the mean across channels and patients to have single average
% band power for every ROI.
% 
% Then it compares the noramtive map to the patient population via z
% scoring and selecting the maximum z-score across the 4 band powers.
% 

% 
% Csaba Kozma
% CNNP Lab, Newcastle University
% July 2023
%% Calcualte mean and sd for control subjects
% Select ROIs based on parcellation
%% Choose the approprite ROIs based on parcellation scheme
if parc==2
    ROIs_normative=normative_table.ROI_2;
    ROIs_patient=patient_table.ROI_2;
elseif parc==4
    ROIs_normative=normative_table.ROI_4;
    ROIs_patient=patient_table.ROI_2;
end

%Assign the number of ROIs and number of subjects
nbroi=max(ROIs_normative);
nbUCLsubj=length(UCLsubjID);

% Preallocate mean and std and subject data
ROIdatamean=nan(nbroi,nbands);
ROIdatastd=nan(nbroi,nbands);
ROIdatanbsubj=nan(nbroi,1);



for i=1:nbroi %nuber of ROIs, changes based on parcellation scheme


roilocs=find(ROIs_normative==i); %find the ROIs

uniquesubj=unique(normative_table.SubjName(roilocs));
currsubjdata_mean=nan(numel(uniquesubj),nbands);
currsubjdata_std=nan(numel(uniquesubj),nbands);
for j=1:numel(uniquesubj)
    [~,locs]=ismember(normative_table.SubjName,uniquesubj{j}); %cycle through control subject by subject
    locs=locs .* ROIs_normative==i;
    locs=find(locs);
    
    currentdata=normative_table.rel_bp(locs,:); %select corresponding channels with bandpower 
    currsubjdata_mean(j,:)=mean(currentdata,1); % take mean across channels
    currsubjdata_std(j,:)=std(currentdata,[],1); % take sd across channels
end
ROIdatamean(i,:)=mean(currsubjdata_mean,1); % take mean across subjects
ROIdatastd(i,:)=std(currsubjdata_mean,[],1);% take sd across subjects
ROIdatanbsubj(i,1)=numel(uniquesubj); %subject number
end
clear currentdata currsubjdata_mean currsubjdata_std currentdatapsd currsubjdatapsd_mean currsubjdatapsd_std i j locs uniquesubj roilocs
disp(['Mean ',num2str(mean(ROIdatanbsubj)),' patients per ROI', ', (std ',num2str(std(ROIdatanbsubj)),').'])

%% Calculate mean and sd for patient population
UCLROIdatamean=nan(nbroi,nbands,nbUCLsubj);
UCLROIdatastd=nan(nbroi,nbands,nbUCLsubj);
UCLROIcontainsresected=ones(nbroi,nbUCLsubj).*-1;

for i=1:nbroi
for j=1:nbUCLsubj

    uclroilocs=find(ROIs_patient==i & (ismember(patient_table.SubjName,str2double(UCLsubjID{j}))));
    if numel(uclroilocs)>0

        uclcurrentdata=patient_table.rel_bp(uclroilocs,:);
        UCLROIdatamean(i,:,j)=mean(uclcurrentdata,1)';% take mean across patients
        UCLROIdatastd(i,:,j)=std(uclcurrentdata,[],1);% take sd across patients
        
        UCLROIcontainsresected(i,j)=mean(patient_table.isResected(uclroilocs)); %select resected ROIs
    end
end
end
clear i j uclroilocs uclcurrentdata uclcurrentdatapsd

%% Calculate abnromality
UCLROIdata_abnormality=nan(size(UCLROIdatamean));

for i=1:nbUCLsubj
    resected=find(UCLROIcontainsresected(:,i)>resectedThresh); %separate resected ROIs
    spared=find(UCLROIcontainsresected(:,i)<=resectedThresh & UCLROIcontainsresected(:,i)~=-1); %separate spared ROIs

    z_UCLH_resected=nan(numel(resected),nbands);
    z_UCLH_spared=nan(numel(spared),nbands);
    for j=1:nbands %for each band
        data1muclh=(UCLROIdatamean(resected,j,i));% mean UCLH BP
        data1suclh=(UCLROIdatastd(resected,j,i));% SD UCLH BP
        data1mnorm=(ROIdatamean(resected,j));% mean Normative BP
        data1snorm=(ROIdatastd(resected,j));% SD Normative BP
        z_UCLH_resected(:,j)=getES(data1muclh,data1suclh,data1mnorm,data1snorm,1);
        
        UCLROIdata_abnormality(resected,j,i)=z_UCLH_resected(:,j);
        
        data2mucl=(UCLROIdatamean(spared,j,i));% mean UCLH BP
        data2sucl=(UCLROIdatastd(spared,j,i));% SD UCLH BP
        data2mnorm=(ROIdatamean(spared,j));% mean Normative BP
        data2snorm=(ROIdatastd(spared,j));% SD Normative BP
        z_UCLH_spared(:,j)=getES(data2mucl,data2sucl,data2mnorm,data2snorm,1); %difference between resected and spared
        
        
        UCLROIdata_abnormality(spared,j,i)=z_UCLH_spared(:,j);
        
        UCLROIdata_abnormality(UCLROIdata_abnormality==inf)=nan;
        clear data1mucl data1sucl data1mnorm data1snorm data2mucl data2sucl data2mnorm data2snorm        
   
    end 
    UCLROIdata_abnormality_tmp=reshape(max(abs(UCLROIdata_abnormality),[],2),nbroi,nbUCLsubj);    
    UCLROIdata_abnormality_final(:,i)=UCLROIdata_abnormality_tmp(:,i);
    
end
end