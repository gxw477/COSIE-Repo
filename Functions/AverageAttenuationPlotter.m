function [f1 ,f2] = AverageAttenuationPlotter(data,plotData,unsegData,cmap)


    nEML = size(plotData,1);

    f1 = figure 
    plot(data{1}.allDepths,-1*data{1}.attSubCut + 10.*log10(unsegData(:,1)),'k-o','MarkerFaceColor','k')
    hold on 
    for i = 1:nEML
        plot(data{1}.allDepths,-1*data{1}.attSubCut + 10.*log10(plotData(i,:,1)),'-o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:))
    end
    plot(data{1}.allDepths,(data{1}.allDepths-data{1}.liverStart*1e3).*0.1*data{1}.liverAtt,'k-.')
    xlabel('Depth (mm)')
    ylabel('Attenuation (dB)')
    set(gca,'FontSize',20)
    xlim([10 60])
    
    
    
    f2 = figure
    plot(data{1}.allDepths,10.*log10(unsegData(:,2)),'k-o','MarkerFaceColor','k')
    hold on 
    for i = 1:nEML
        plot(data{1}.allDepths, 10.*log10(plotData(i,:,2)),'-o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:))
    end
    
    xlabel('Depth (mm)')
    ylabel('RMSE (dB)')
    xlim([10 60])
    set(gca,'FontSize',20)
    colormap([cmap]);
    
    colormap([0,0,0;cmap])
    cB = colorbar;
    set(cB,'XTick',(1:2:(2*(nEML+1)))./(2*(nEML+1)))
    set(cB,'TickLabels',0:nEML)



end
