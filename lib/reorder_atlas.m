function [atlas_tbl]=reorder_atlas()

%% Reorder Atlas

% atlas files
file_atlas = 'atlas/ATLAS.mat';
file_retains = 'atlas/retains.mat';
file_atlas_info = 'atlas/atlasinfo.mat';


%% load atlas
load(file_retains)
load(file_atlas)
load(file_atlas_info)

%% combine info and format

scale = {'36','60','125','250'}';

parc_name = strcat('scale',scale);
n_parc = length(parc_name);

name = cell(n_parc,1);
xyz = cell(n_parc,1);
vol = cell(n_parc,1);
dists = cell(n_parc,1);
retain = cell(n_parc,1);
n_roi = zeros(n_parc,1);

for i=1:n_parc
    name{i} = eval(['ATLAS.name' scale{i}]);
    n_roi(i) = length(name{i});
    xyz{i} = eval(['xyz' scale{i}]);
    vol{i} = eval(['vol' scale{i}]);
    dists{i} = eval(['dists' scale{i}]);
    retain{i} = eval(['scale' scale{i} 'retain']);
end


%% table version for MATLAB
atlas_tbl = table(parc_name,n_roi,name,xyz,vol,dists,retain);

%% Getting the retained ROIs
for i=1:size(atlas_tbl.retain,1)

atlas_tbl.name{i}=atlas_tbl.name{i}(atlas_tbl.retain{i}==1,:);
atlas_tbl.xyz{i}=atlas_tbl.xyz{i}(atlas_tbl.retain{i}==1,:);
atlas_tbl.vol{i}=atlas_tbl.vol{i}(atlas_tbl.retain{i}==1,:);
atlas_tbl.dists{i}=atlas_tbl.dists{i}(atlas_tbl.retain{i}==1,:);

end
clear i
    

end