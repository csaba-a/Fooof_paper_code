function save_figs(patient,type_power_spectrum,figdir)
%% Export figs
patient= char(patient);
if ~contains(patient, 'Group')
    Fig_labels = {['Brain_abnormality_top_',type_power_spectrum, '_',patient], ['DRS_',type_power_spectrum,'_cmp_',patient]};

else
    Fig_labels = {['Beeswarm_plot_',type_power_spectrum, '_cmp_',patient], ['AUC_',type_power_spectrum,'_cmp_',patient]};
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    fig_idx = get(FigHandle, 'Number');
    FigName   =  [Fig_labels{fig_idx}];

    set(0, 'CurrentFigure', FigHandle);

    saveas(FigHandle, fullfile(figdir,strcat(FigName, '.pdf'))); % specify the full path
end
close all

end