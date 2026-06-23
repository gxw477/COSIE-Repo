function [powerUnsegMean2,powerUnsegStd2,powerSegMean2,powerSegStd2] = doubleSegmentationAveragePlot(powerUnseg,powerSegTest,jitter,xPlot,mubsTH,mubsStd,meanIdx,depthAverageIdxs,imageAverageIdxs,cmap)
    
    powerUnseg = powerUnseg(imageAverageIdxs,depthAverageIdxs,:);
    powerSegTest = powerSegTest(imageAverageIdxs,depthAverageIdxs,:,:); 

    fillColor = [140, 222, 162]./256;

    pos = mubsStd(2);
    neg = mubsStd(1);

    nNormUnseg = sum(~isnan(powerUnseg(:,:,1)),meanIdx);
    powerUnsegMean2 = squeeze(mean(powerUnseg(:,:,1),meanIdx,'omitnan'));
    powerUnsegStd2  = squeeze(1./nNormUnseg.*sqrt(sum(powerUnseg(:,:,2).^2,meanIdx,'omitnan')));
    
    nNorm = sum(~isnan(powerSegTest(:,:,1)),meanIdx);
    powerSegMean2 = squeeze(mean(powerSegTest(:,:,:,1),meanIdx,'omitnan'));
    powerSegStd2 = squeeze(1./nNorm.*sqrt(sum(powerSegTest(:,:,:,2).^2,meanIdx,'omitnan')));
    
    powerSegStd2(powerSegStd2==inf) = .1;
    %jitter = linspace(-.3,.3,nEML+1);
    
    figure 
    xShade = [0 0 65 65 ];
    yShade = [ neg pos pos neg ];
    area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',neg,'ShowBaseLine','off');
    hold on 
    plot([xShade(1) xShade(end)],mubsTH.*[1,1],'-.','Color',[255 171 0]./256,'LineWidth',2)
    plot([xShade(1) xShade(end)],pos.*[1 1],'-.','LineWidth',2,'Color',fillColor)
    plot([xShade(1) xShade(end)],neg.*[1 1],'-.','LineWidth',2,'Color', fillColor)
    plot([xShade(1) xShade(1)],[neg pos],'-.','LineWidth',2,'Color',fillColor)
    plot([xShade(end) xShade(end)],[neg pos],'-.','LineWidth',2,'Color',fillColor)
    
    
    errorbar((xPlot)+jitter(1),powerUnsegMean2,powerUnsegStd2,'k.','MarkerSize',12)
    hold on 
    for i = 1:size(powerSegMean2,2)
        errorbar((xPlot) + jitter(i+1) ,powerSegMean2(:,i),powerSegStd2(:,i),'.','Color',cmap(i,:),'MarkerSize',12)
    end
    
    set(gca,'YScale','log')
    set(gca,'FontSize',20)
    ylim([1e-2 1e1])
        
    if meanIdx == 1 
        xlabel('Depth (mm)')
        xlim([20 60])
        set(gca,'XTick',[20:10:60])
    elseif meanIdx == 2 
        xlabel('Image')
        xlim([0 8])
        set(gca,'XTick',[1:7])
    end

    ylabel('BSC (m^{-1}sr^{-1})')
   

end
