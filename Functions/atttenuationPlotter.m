function [] = atttenuationPlotter(allDepths,unsegData,segArray,liverStart,liverAtt,attSubCut)
    
    
    nEML = size(segArray,1);
    cmap = [
            0.1216 0.4667 0.7059;  % Strong Blue
            0.8510 0.3725 0.0078;  % Vivid Orange
            0.5569 0.2667 0.6784; 
            0.7922 0.1176 0.5686;  % Magenta
            0.4980 0.4980 0.4980;  % Dark Gray
        ];

    figure
    plot(allDepths,-attSubCut + 10.*log10(unsegData(:,1)),'k-o','MarkerFaceColor','k')
    hold on 
    for i = 1:nEML
        plot(allDepths,-attSubCut + 10.*log10(segArray(i,:,1)),'-o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:))
    end
    
    plot(allDepths,(allDepths-liverStart*1e3).*0.1*liverAtt,'k-.')
    xlabel('Depth (mm)')
    ylabel('Attenuation (dB)')
    set(gca,'FontSize',20)
    xlim([10 60])
    
    
    
    figure
    plot(allDepths,10.*log10(unsegData(:,2)),'k-o','MarkerFaceColor','k')
    hold on 
    for i = 1:nEML
        plot(allDepths,10.*log10(segArray(i,:,2)),'-o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:))
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
