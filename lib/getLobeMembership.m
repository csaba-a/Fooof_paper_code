function [lobemembership,lobes]=getLobeMembership(T,parc)
lobes={'l.subcortical','l.temporal','l.frontal','l.parietal','l.occipital',...
    'r.subcortical','r.temporal','r.frontal','r.parietal','r.occipital'};

lobemembership=nan(size(T.name{parc},1),1);

subcort=1:14;
lobemembership(1:7)=1;
% Left
% temporal 2
tmp=[find(contains(T.name{parc},'l.superiortemporal'))];
tmp=[tmp;find(contains(T.name{parc},'l.middletemporal'))];
tmp=[tmp;find(contains(T.name{parc},'l.inferiortemporal'))];
tmp=[tmp;find(contains(T.name{parc},'l.bankssts'))];
tmp=[tmp;find(contains(T.name{parc},'l.fusiform'))];
tmp=[tmp;find(contains(T.name{parc},'l.transversetemporal'))];
tmp=[tmp;find(contains(T.name{parc},'l.entorhinal'))];
tmp=[tmp;find(contains(T.name{parc},'l.temporalpole'))];
tmp=[tmp;find(contains(T.name{parc},'l.parahippocampal'))];
tmp=[tmp;find(contains(T.name{parc},'l.insula'))];
lobemembership(tmp)=2;

% Left
% frontal 3
tmp=[find(contains(T.name{parc},'l.superiorfrontal'))];
tmp=[tmp;find(contains(T.name{parc},'l.rostralmiddlefrontal'))];
tmp=[tmp;find(contains(T.name{parc},'l.caudalmiddlefrontal'))];
tmp=[tmp;find(contains(T.name{parc},'l.parsopercularis'))];
tmp=[tmp;find(contains(T.name{parc},'l.parsorbitalis'))];
tmp=[tmp;find(contains(T.name{parc},'l.parstriangularis'))];
tmp=[tmp;find(contains(T.name{parc},'l.lateralorbitofrontal'))];
tmp=[tmp;find(contains(T.name{parc},'l.medialorbitofrontal'))];
tmp=[tmp;find(contains(T.name{parc},'l.precentral'))];
tmp=[tmp;find(contains(T.name{parc},'l.paracentral'))];
tmp=[tmp;find(contains(T.name{parc},'l.frontalpole'))];
tmp=[tmp;find(contains(T.name{parc},'l.caudalanteriorcingulate'))];
tmp=[tmp;find(contains(T.name{parc},'l.rostralanteriorcingulate'))];
lobemembership(tmp)=3;

% Left
% parietal 4
tmp=[find(contains(T.name{parc},'l.superiorparietal'))];
tmp=[tmp;find(contains(T.name{parc},'l.inferiorparietal'))];
tmp=[tmp;find(contains(T.name{parc},'l.supramarginal'))];
tmp=[tmp;find(contains(T.name{parc},'l.postcentral'))];
tmp=[tmp;find(contains(T.name{parc},'l.precuneus'))];
tmp=[tmp;find(contains(T.name{parc},'l.isthmuscingulate'))];
tmp=[tmp;find(contains(T.name{parc},'l.posteriorcingulate'))];
lobemembership(tmp)=4;

% Left
% occipital 5
tmp=[find(contains(T.name{parc},'l.lateraloccipital'))];
tmp=[tmp;find(contains(T.name{parc},'l.lingual'))];
tmp=[tmp;find(contains(T.name{parc},'l.cuneus'))];
tmp=[tmp;find(contains(T.name{parc},'l.pericalcarine'))];
lobemembership(tmp)=5;


lobemembership(8:14)=6;
% r
% temporal 7
tmp=[find(contains(T.name{parc},'r.superiortemporal'))];
tmp=[tmp;find(contains(T.name{parc},'r.middletemporal'))];
tmp=[tmp;find(contains(T.name{parc},'r.inferiortemporal'))];
tmp=[tmp;find(contains(T.name{parc},'r.bankssts'))];
tmp=[tmp;find(contains(T.name{parc},'r.fusiform'))];
tmp=[tmp;find(contains(T.name{parc},'r.transversetemporal'))];
tmp=[tmp;find(contains(T.name{parc},'r.entorhinal'))];
tmp=[tmp;find(contains(T.name{parc},'r.temporalpole'))];
tmp=[tmp;find(contains(T.name{parc},'r.parahippocampal'))];
tmp=[tmp;find(contains(T.name{parc},'r.insula'))];
lobemembership(tmp)=7;

% r
% frontal 8
tmp=[find(contains(T.name{parc},'r.superiorfrontal'))];
tmp=[tmp;find(contains(T.name{parc},'r.rostralmiddlefrontal'))];
tmp=[tmp;find(contains(T.name{parc},'r.caudalmiddlefrontal'))];
tmp=[tmp;find(contains(T.name{parc},'r.parsopercularis'))];
tmp=[tmp;find(contains(T.name{parc},'r.parsorbitalis'))];
tmp=[tmp;find(contains(T.name{parc},'r.parstriangularis'))];
tmp=[tmp;find(contains(T.name{parc},'r.lateralorbitofrontal'))];
tmp=[tmp;find(contains(T.name{parc},'r.medialorbitofrontal'))];
tmp=[tmp;find(contains(T.name{parc},'r.precentral'))];
tmp=[tmp;find(contains(T.name{parc},'r.paracentral'))];
tmp=[tmp;find(contains(T.name{parc},'r.frontalpole'))];
tmp=[tmp;find(contains(T.name{parc},'r.caudalanteriorcingulate'))];
tmp=[tmp;find(contains(T.name{parc},'r.rostralanteriorcingulate'))];
lobemembership(tmp)=8;
% r
% parietal 9
tmp=[find(contains(T.name{parc},'r.superiorparietal'))];
tmp=[tmp;find(contains(T.name{parc},'r.inferiorparietal'))];
tmp=[tmp;find(contains(T.name{parc},'r.supramarginal'))];
tmp=[tmp;find(contains(T.name{parc},'r.postcentral'))];
tmp=[tmp;find(contains(T.name{parc},'r.precuneus'))];
tmp=[tmp;find(contains(T.name{parc},'r.isthmuscingulate'))];
tmp=[tmp;find(contains(T.name{parc},'r.posteriorcingulate'))];
lobemembership(tmp)=9;

% r
% occipital 10
tmp=[find(contains(T.name{parc},'r.lateraloccipital'))];
tmp=[tmp;find(contains(T.name{parc},'r.lingual'))];
tmp=[tmp;find(contains(T.name{parc},'r.cuneus'))];
tmp=[tmp;find(contains(T.name{parc},'r.pericalcarine'))];
lobemembership(tmp)=10;

%%
end