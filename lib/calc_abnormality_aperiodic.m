function [UCLROIdata_abnormality_aperiodic,UCLROIcontainsresected]=calc_abnormality_aperiodic(normative_table, patient_table,UCLsubjID,resectedThresh, parc)

%% NORMATIVE code averages within a patient first, then across patients
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

ROIdatamean=nan(nbroi,1);
ROIdatastd=nan(nbroi,1);

for i=1:nbroi


    roilocs=find(ROIs_normative(:)==i);

    currentdata=normative_table.aperiodic_cmps_2(roilocs,:);

    ROIdatamean(i,:)=mean(currentdata,1);
    ROIdatastd(i,:)=std(currentdata,[],1);
end
clear currentdata currsubjdata_mean currsubjdata_std currentdatapsd currsubjdatapsd_mean currsubjdatapsd_std i j locs uniquesubj roilocs
%%
UCLROIdatamean=nan(nbroi,nbUCLsubj);
UCLROIdatastd=nan(nbroi,nbUCLsubj);
UCLROIcontainsresected=ones(nbroi,nbUCLsubj).*-1;

for i=1:nbroi
    for j=1:nbUCLsubj

        uclroilocs=find(ROIs_patient(:)==i & (strcmp(string(patient_table.SubjName(:)),UCLsubjID{j})));
        if numel(uclroilocs)>0
            %disp(string(uclroilocs))
            uclcurrentdata=patient_table.aperiodic_cmps_2(uclroilocs);
            UCLROIdatamean(i,j)=mean(uclcurrentdata,1)';
            UCLROIdatastd(i,j)=std(uclcurrentdata,[],1);

            UCLROIcontainsresected(i,j)=mean(patient_table.isResected(uclroilocs));
        end
    end
end
clear i j uclroilocs uclcurrentdata uclcurrentdatapsd

%% Calculate abnromality
UCLROIdata_abnormality_aperiodic=nan(size(UCLROIdatamean));


for i=1:nbUCLsubj
    resected=find(UCLROIcontainsresected(:,i)>resectedThresh);
    spared=find(UCLROIcontainsresected(:,i)<=resectedThresh & UCLROIcontainsresected(:,i)~=-1);

    z_UCLH_resected=nan(numel(resected),1);
    z_UCLH_spared=nan(numel(spared),1);
    %for each band
    data1muclh=(UCLROIdatamean(resected,i));% mean UCLH BP
    data1suclh=(UCLROIdatastd(resected,i));% SD UCLH BP
    data1mnorm=(ROIdatamean(resected));% mean Normative BP
    data1snorm=(ROIdatastd(resected));% SD Normative BP
    z_UCLH_resected(:)=getES(data1muclh,data1suclh,data1mnorm,data1snorm,1);

    UCLROIdata_abnormality_aperiodic(resected,i)=abs(z_UCLH_resected(:));

    data2mucl=(UCLROIdatamean(spared,i));% mean UCLH BP
    data2sucl=(UCLROIdatastd(spared,i));% SD UCLH BP
    data2mnorm=(ROIdatamean(spared));% mean Normative BP
    data2snorm=(ROIdatastd(spared));% SD Normative BP
    z_UCLH_spared(:)=getES(data2mucl,data2sucl,data2mnorm,data2snorm,1);


    UCLROIdata_abnormality_aperiodic(spared,i)=abs(z_UCLH_spared(:));

    UCLROIdata_abnormality_aperiodic(UCLROIdata_abnormality_aperiodic==inf)=nan;
    clear data1mucl data1sucl data1mnorm data1snorm data2mucl data2sucl data2mnorm data2snorm


end
end