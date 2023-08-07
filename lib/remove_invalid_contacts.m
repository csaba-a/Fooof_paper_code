%% Remove invalid data and disregard >30Hz data
% remove invalid contacts (e.g. in white matter etc)
toremove=sum(isnan(MasterChannelTable.ROIID_inAtlas),2)>0;

disp(strcat('Removing  ',num2str(sum(toremove)),' contacts (',num2str(100*(sum(toremove)/numel(toremove))),'%).'));
MasterChannelTable(toremove==1,:)=[];
clear toremove

%only need up to 30Hz
MasterChannelTable.pxx_n(:,60:end)=[]; %0.5 overlap on the frequency dimension means 1-60 datapoints is 1-30Hz


%% Define boolian for separating datasets

UCLH_bool = MasterChannelTable.DataSource=="UCLH";
RAM_bool = MasterChannelTable.DataSource=="RAM";
