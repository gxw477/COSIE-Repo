function [] = bscSegPlotterDual(powerSegCOH,powerSEGSNR,convFactor,mubsTH,mubsTH_STD,logSwitch)


    %% Coherence Section
    
    [~,uIDs,~] = unique(powerSegCOH(4,:));

    powerSegCOH = powerSegCOH(:,uIDs);
    
    fillColor = [140, 222, 162]./256;

    bscAllMean_COH = powerSegCOH(1,:).*convFactor;
    bscAllStd_COH = powerSegCOH(2,:).*convFactor;

    bscAllStdNeg_COH = bscAllStd_COH;
    
    lTHbool_COH= bscAllMean_COH-bscAllStd_COH < 0;
    bscAllStdNeg_COH(lTHbool_COH) = bscAllMean_COH(lTHbool_COH);

    %% SNR Section
    

    
    [~,uIDs,~] = unique(powerSEGSNR(4,:));

    powerSEGSNR = powerSEGSNR(:,uIDs);
    
    %fillColor = [140, 222, 162]./256;

    bscAllMean_SNR = powerSEGSNR(1,:).*convFactor;
    bscAllStd_SNR = powerSEGSNR(2,:).*convFactor;

    bscAllStdNeg_SNR = bscAllStd_SNR;
    
    lTHbool_SNR = bscAllMean_SNR-bscAllStd_SNR < 0;
    bscAllStdNeg_SNR(lTHbool_SNR) = bscAllMean_SNR(lTHbool_SNR);

    
    
    if length(mubsTH_STD) > 1 
        neg = mubsTH_STD(1);
        pos = mubsTH_STD(2);
    else
        error('Dickhead')
    end


    %%

    figure 
    xShade = [-15 -15 110 110 ];
    yShade = [ neg pos pos neg ];
    area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',neg,'ShowBaseLine','off');
    
    
    area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',neg,'ShowBaseLine','off');
    hold on 
    plot([-15 110],mubsTH.*[1,1],'-.','Color',[255 171 0]./256,'LineWidth',2)
    plot([-15 110],pos.*[1 1],'-.','LineWidth',2,'Color',fillColor)
    plot([-15 110],neg.*[1 1],'-.','LineWidth',2,'Color', fillColor)
    plot([-15 -15],[neg pos],'-.','LineWidth',2,'Color',fillColor)
    plot([110 110],[neg pos],'-.','LineWidth',2,'Color',fillColor)
    

    

    errorbar(powerSegCOH(4,1),bscAllMean_COH(1),bscAllStdNeg_COH(1),bscAllStd_COH(1),'r.','MarkerSize',10) 
    errorbar(powerSegCOH(4,2:2:end),bscAllMean_COH(2:2:end),bscAllStdNeg_COH(2:2:end),bscAllStd_COH(2:2:end),'k-o','MarkerSize',5,'MarkerFaceColor','k') 
    errorbar(powerSEGSNR(4,2:2:end),bscAllMean_SNR(2:2:end),bscAllStdNeg_SNR(2:2:end),bscAllStd_SNR(2:2:end),'-o','Color',[0 2 255]./256,'MarkerFaceColor','blue') 

    %errorbar(-10,mubsTH,mubsTH_STD,'r.')
    if logSwitch
        set(gca,'YScale','log')
        yticks(10.^([-4:1]))
        set(gca,'YLim',[1e-4 10])
    end

    %yticks(0.05.*2.^[-1:6])

    xticks([0 : 25: 100])
    xlabel('Segmentation (%)')
    ylabel('BSC (m^{-1}sr^{-1})')
    set(gca,'FontSize',20)
    %ylim padded
    xlim([-14 104])


    
end

