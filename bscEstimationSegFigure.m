function [] = bscEstimationSegFigure(bscSpeckleBf,bscSpeckleSTD,bscEstimate)

    
    figure
    tL = tiledlayout(1,1);
    ax1 = axes(tL);
    
    errorbar(ax1,-10, bscSpeckleBf,bscSpeckleSTD,'r.')
    hold on 
    plot(ax1,[-10 bscEstimate(end,3)+10],bscSpeckleBf.*[1 1],'r-.')
    errorbar(ax1,bscEstimate(1,3),bscEstimate(1,1),bscEstimate(1,2),'k.')
    errorbar(ax1,bscEstimate(2:end,3),bscEstimate(2:end,1),bscEstimate(2:end,2),'k.')
    xlabel(ax1,'Segmentation (%)')
    ylabel(ax1,'BSC (m^{-1}sr^{-1})')
    xticks(ax1,[0:10:80])
    xlim(ax1,[-20 bscEstimate(end,3)+10])
    yLim2 = max([bscEstimate(2:end,1)+bscEstimate(2:end,2); bscSpeckleBf+bscSpeckleSTD] )*1.05;
    ylim(ax1,[0 yLim2])
    
    ax2 = axes(tL);
    plot(ax2,bscEstimate(2:end,3),bscEstimate(2:end,1),'k.')
    ylim(ax2,[0 yLim2])
    xlim(ax2,[-20 bscEstimate(end,3)+10])
    xticks(ax2,unique(bscEstimate(:,3)))
    yticks(ax2,[])
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'left';
    set(ax2,'YColor','k')
    set(ax2,'XColor','k')
    
    xticklabels(ax2,num2cell(round(bscEstimate(:,4))))
    set(ax2,'FontSize',8)
    ax2.Color = 'none';
    ax1.Box = 'off';
    ax2.Box = 'off';

end