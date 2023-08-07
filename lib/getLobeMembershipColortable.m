function [label,lobes]=getLobeMembershipColortable(colortable,label)
lobes={'temporal','frontal','parietal','occipital'};

colortable.lobe=nan(size(colortable.struct_names));

% temporal 1
tmp=[find(contains(colortable.struct_names,'superiortemporal'))];
tmp=[tmp;find(contains(colortable.struct_names,'middletemporal'))];
tmp=[tmp;find(contains(colortable.struct_names,'inferiortemporal'))];
tmp=[tmp;find(contains(colortable.struct_names,'bankssts'))];
tmp=[tmp;find(contains(colortable.struct_names,'fusiform'))];
tmp=[tmp;find(contains(colortable.struct_names,'transversetemporal'))];
tmp=[tmp;find(contains(colortable.struct_names,'entorhinal'))];
tmp=[tmp;find(contains(colortable.struct_names,'temporalpole'))];
tmp=[tmp;find(contains(colortable.struct_names,'parahippocampal'))];

tmp=[tmp;find(contains(colortable.struct_names,'insula'))];

colortable.lobe(tmp)=1;

% frontal 2
tmp=[find(contains(colortable.struct_names,'superiorfrontal'))];
tmp=[tmp;find(contains(colortable.struct_names,'rostralmiddlefrontal'))];
tmp=[tmp;find(contains(colortable.struct_names,'caudalmiddlefrontal'))];
tmp=[tmp;find(contains(colortable.struct_names,'parsopercularis'))];
tmp=[tmp;find(contains(colortable.struct_names,'parsorbitalis'))];
tmp=[tmp;find(contains(colortable.struct_names,'parstriangularis'))];
tmp=[tmp;find(contains(colortable.struct_names,'lateralorbitofrontal'))];
tmp=[tmp;find(contains(colortable.struct_names,'medialorbitofrontal'))];
tmp=[tmp;find(contains(colortable.struct_names,'precentral'))];
tmp=[tmp;find(contains(colortable.struct_names,'paracentral'))];
tmp=[tmp;find(contains(colortable.struct_names,'frontalpole'))];
tmp=[tmp;find(contains(colortable.struct_names,'caudalanteriorcingulate'))];
tmp=[tmp;find(contains(colortable.struct_names,'rostralanteriorcingulate'))];
colortable.lobe(tmp)=2;

% parietal 3
tmp=[find(contains(colortable.struct_names,'superiorparietal'))];
tmp=[tmp;find(contains(colortable.struct_names,'inferiorparietal'))];
tmp=[tmp;find(contains(colortable.struct_names,'supramarginal'))];
tmp=[tmp;find(contains(colortable.struct_names,'postcentral'))];
tmp=[tmp;find(contains(colortable.struct_names,'precuneus'))];
tmp=[tmp;find(contains(colortable.struct_names,'isthmuscingulate'))];
tmp=[tmp;find(contains(colortable.struct_names,'posteriorcingulate'))];
colortable.lobe(tmp)=3;

% occipital 4
tmp=[find(contains(colortable.struct_names,'lateraloccipital'))];
tmp=[tmp;find(contains(colortable.struct_names,'lingual'))];
tmp=[tmp;find(contains(colortable.struct_names,'cuneus'))];
tmp=[tmp;find(contains(colortable.struct_names,'pericalcarine'))];
colortable.lobe(tmp)=4;

% % left cingulate 5
% tmp=[find(contains(colortable.struct_names,'isthmuscingulate'))];
% tmp=[tmp;find(contains(colortable.struct_names,'posteriorcingulate'))];
% tmp=[tmp;find(contains(colortable.struct_names,'caudalanteriorcingulate'))];
% tmp=[tmp;find(contains(colortable.struct_names,'rostralanteriorcingulate'))];
% colortable.lobe(tmp)=5;

% left insula
% tmp=[find(contains(colortable.struct_names,'insula'))];
% colortable.lobe(tmp)=6;

for i=1:numel(label)
    label(i)=colortable.lobe(find(colortable.table(:,5)==label(i)));
end



%%
end