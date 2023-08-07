function makeFigNice()
%     if hFig==[]
        hFig=gcf;
%     end
    ax=get(hFig,'CurrentAxes');
    set(ax,'FontSize',18);
    set(hFig,'Color',[1 1 1])
%     set(findobj('Type','line'),'Color','k')
    
end