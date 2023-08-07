function [aucdrs]=calc_auc(UCLROIdata_abnormality,UCLROIcontainsresected,UCLsubjID,resectedThresh)
%% CALC_AUC
% Calculates auc (DRS) values for all subjects that are imported (individual or group level)
%
%%
%Preallocate auc varaible
aucdrs=nan(length(UCLsubjID),1);

for i=1:length(UCLsubjID)
    resected=find(UCLROIcontainsresected(:,i)>resectedThresh); %separate resected ROIs
    spared=find(UCLROIcontainsresected(:,i)<=resectedThresh & UCLROIcontainsresected(:,i)~=-1); %separate spared ROIs

    
    try % calculate AUC
        [~,~,~,aucdrs(i,1)]=perfcurve([zeros(size(UCLROIdata_abnormality(resected,i)));ones(size(UCLROIdata_abnormality(spared,i)))],([UCLROIdata_abnormality(resected,i);UCLROIdata_abnormality(spared,i)]),1);
    catch
    end
end
