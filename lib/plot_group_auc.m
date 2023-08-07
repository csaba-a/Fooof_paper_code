function plot_group_auc(aucdrs, metadata)

% calculate AUC
clc
disp('Using normative map')
[~,p]=ttest2(aucdrs(find(metadata.ILAE1<3)),aucdrs(find(metadata.ILAE1>2)),'tail','left')
[~,~,~,auc]=perfcurve([zeros(size(find(metadata.ILAE1<3)));ones(size(find(metadata.ILAE1>2)))],[aucdrs(find(metadata.ILAE1<3));aucdrs(find(metadata.ILAE1>2))],1)

% plot beeswarm plot
figure
set(gcf,'renderer','painters');
violin(padWithNans(aucdrs(find(metadata.ILAE1<3)),aucdrs(find(metadata.ILAE1>2))),'facecolor',[1 1 1]);
UnivarScatter(padWithNans(aucdrs(find(metadata.ILAE1<3)),nan(1,2)),'Label',{'ILAE1,2','ILAE3+'},'MarkerFaceColor','blue','MarkerEdgeColor','none');

hold on
UnivarScatter(padWithNans(nan(1,2),aucdrs(find(metadata.ILAE1>2))),'Label',{'ILAE1,2','ILAE3+'},'MarkerFaceColor','red','MarkerEdgeColor','none');
plot([0.5 2.5],[0.5 0.5],'k:');
ylabel('D_R_S')
box off
ylim([0 1])
set(gcf,'position', [10 10 300 300])
title(['p=',num2str(p,2)],[' AUC= ', num2str(auc,2)])


% plot AUC
[X,Y,T,auc,optrocpt,subY,subYnames]=perfcurve([zeros(size(find(metadata.ILAE1<3)));ones(size(find(metadata.ILAE1>2)))],[aucdrs(find(metadata.ILAE1<3));aucdrs(find(metadata.ILAE1>2))],1);
figure
set(gcf,'renderer','painters');
plot(X,Y,'k','linewidth',2)
hold on
plot([0 1],[0 1],'r','linewidth',2)
axis square
makeFigNice
scatter(optrocpt(1),optrocpt(2),100,'filled','k')
title(['p= ',num2str(p),' AUC= ', num2str(auc)])

end
