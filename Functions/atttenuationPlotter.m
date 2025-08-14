function [] = atttenuationPlotter(allDepths,unsegData,segArray,liverStart,liverAtt,attSubCut)
    
    
    nEMLidx = size(segArray,1);
    greenValues = round(linspace(255,0,nEMLidx));
    cmap = [255.*ones(nEMLidx,1),greenValues',zeros(nEMLidx,1)]./255;
  
    figure 
   
    plot(allDepths,attSubCut + 10.*log10(unsegData(:,1)),'k-o','MarkerFaceColor','k')
    hold on 
    for i = 1:nEMLidx
        plot(allDepths,attSubCut + 10.*log10(segArray(i,:,1)),'-o','Color',cmap(i,:))
    end
    
    plot(allDepths,attSubCut + (allDepths-liverStart*1e3).*1e-3*liverAtt,'k-.')
    xlabel('Depth (mm)')
    ylabel('Attenuation (dB)')
    set(gca,'FontSize',20)
    xlim padded
    
    
    
    figure
    plot(allDepths,10.*log10(unsegData(:,2)),'k-o','MarkerFaceColor','k')
    hold on 
    for i = 1:nEMLidx
        plot(allDepths,10.*log10(segArray(i,:,2)),'-o','Color',cmap(i,:))
    end
    
    xlabel('Depth (mm)')
    ylabel('RMSE (dB)')
    xlim padded
    set(gca,'FontSize',20)
    colormap(cmap)
    cB = colorbar;
    set(cB,'XTick',0.05:0.1:1)
    set(cB,'XTickLabel',1:10)
    %saveas(gcf,[imgDir,'\AttCorr.fig'])
    %exportgraphics(gcf,[imgDir,'\AttCorr.pdf'])

end
