function [MasterChannelTable]=reorder_ieeg(MasterChannelTable)

% Get rid of nan rows
MasterChannelTable.ROIID_inAtlas=double(MasterChannelTable.ROIID_inAtlas);
MasterChannelTable.ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas==0)=nan;

%% directories/files/fields

% atlas files
file_atlas = 'atlas/ATLAS.mat';
file_retains = 'atlas/retains.mat';
file_atlas_info = 'atlas/atlasinfo.mat';


%% load atlas
load(file_retains)
load(file_atlas)
load(file_atlas_info)
%% Prune ROIs from data table
count=1;
ROIID_inAtlas=nan(size(MasterChannelTable.ROIID_inAtlas,1),4);
for i=1:size(scale36retain)
    if scale36retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,1)==i,1)=count;
        count=count+1;
    end
end
count=1;
for i=1:size(scale60retain)
    if scale60retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,2)==i,2)=count;
        count=count+1;
    end
end
count=1;
for i=1:size(scale125retain)
    if scale125retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,3)==i,3)=count;
        count=count+1;
    end
end
count=1;
for i=1:size(scale250retain)
    if scale250retain(i)==1
        ROIID_inAtlas(MasterChannelTable.ROIID_inAtlas(:,4)==i,4)=count;
        count=count+1;
    end
end
MasterChannelTable.ROIID_inAtlas=ROIID_inAtlas;
clear i ROIID_inAtlas count