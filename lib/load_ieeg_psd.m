

%% Preallocate table
DataSource={};
SubjectID={};
ChannelName={};
Channelxyz=[[],[],[]];
Channelside={};
IsResected=[];
pxx=[];
pxx_n=[];
ROIID_inAtlas=[];
MasterChannelTable=table(DataSource,SubjectID,ChannelName,ROIID_inAtlas,Channelxyz,Channelside,IsResected,pxx,pxx_n);
clear DataSource SubjectID ChannelName ROIID_inAtlas Channelxyz Channelside IsResected pxx pxx_n

%% Load UCLH data

load('data/ChannelTable_UCLH.mat')

startlocs=size(MasterChannelTable,1)+1;
endlocs=size(ChannelTable,1)+(startlocs-1);

MasterChannelTable.DataSource(startlocs:endlocs)=cellstr(ChannelTable.DataSource);
MasterChannelTable.SubjectID(startlocs:endlocs)=convertStringsToChars(ChannelTable.SubjectID);
MasterChannelTable.ChannelName(startlocs:endlocs)=cellstr(ChannelTable.ChannelName);
MasterChannelTable.ROIID_inAtlas(startlocs:endlocs,1:4)=ChannelTable.ROIID_inAtlas;

MasterChannelTable.Channelxyz(startlocs:endlocs,1:3)=ChannelTable.Channelxyz;
MasterChannelTable.Channelside(startlocs:endlocs)=cellstr(ChannelTable.Channelside);
MasterChannelTable.IsResected(startlocs:endlocs)=ChannelTable.isResected;
MasterChannelTable.pxx(startlocs:endlocs,1:159)=ChannelTable.pxx;
MasterChannelTable.pxx_n(startlocs:endlocs,1:159)=ChannelTable.pxx_n;

clear ChannelTable startlocs endlocs

%% Load RAM data

load('data/ChannelTable_RAM.mat')

startlocs=size(MasterChannelTable,1)+1;
endlocs=size(ChannelTable,1)+(startlocs-1);
MasterChannelTable.DataSource(startlocs:endlocs)=cellstr(ChannelTable.DataSource);
MasterChannelTable.SubjectID(startlocs:endlocs)=convertStringsToChars(ChannelTable.SubjectID);
MasterChannelTable.ChannelName(startlocs:endlocs)=cellstr(ChannelTable.ChannelName);
MasterChannelTable.ROIID_inAtlas(startlocs:endlocs,1:4)=ChannelTable.ROIID_inAtlas;
MasterChannelTable.pxx_n(startlocs:endlocs,1:159)=ChannelTable.pxx_n;
clear ChannelTable startlocs endlocs

