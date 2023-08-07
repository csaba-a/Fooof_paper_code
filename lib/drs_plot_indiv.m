function []=drs_plot_indiv(resected,spared,maxz)
%% DRS plotting and calcualtion
% calcualtes Drs value, which is the distinguisibility between resected and
% spared region and plots the individual pateint plot
% Csaba Kozma
% CNNP Lab, Newcastle University
% July 2023
%% load colormaps
load('lib/colourmaps_yw.mat'); %this loads cmaps variable

%% Plot
[~,~,~,drs]=perfcurve([zeros(size(resected));ones(size(spared))],maxz([resected;spared]),1);
figure
set(gcf,'position', [10 10 300 300])
set(gcf,'renderer','painters');
hold on
[x,y]=UnivarScatter(padWithNans(maxz(resected),maxz(spared)),'Labels',{'resected','spared'},'MarkerFaceColor',cmaps.c);
x(isnan(x))=[];
y(isnan(y))=[];

c=ceil((y./4)*256);
c(isnan(c))=[];
hold on
scatter(x(:),y(:),ones(size(y(:),1),1)*200,cmaps.z_score_map(c(:),:),'filled')
            
xlim([.5 2.5])
ylim([0 4])
caxis([0 4])
title(strcat('D_R_S=',num2str(drs)))
ylabel('max |z|')
box off
end