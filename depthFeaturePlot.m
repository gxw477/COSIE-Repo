function [] = depthFeaturePlot(cosieStructCOH,cosieStructSNR,EMLidx,lWidth,cohTest,snrTest,sumIdx)

    % depthFeaturePlot(cosieStructCOH,cosieStructSNR,EMLidx,lWidth,cohTest,snrTest)
    % cosieStructCOH : cosie struct containing EML for coherence seg
    % cosieStructSNR : cosie struct containing EML for snr seg
    % EMLidx : self explanatory
    % lWidth : image line width
    % cohTest : coherence values in image 
    % snrTest : snr values for image


    
    figure 
    subplot(2,2,1)
    plot(lWidth.*[1 length(cohTest)],cosieStructCOH.redEML(1,EMLidx).*[1 1],'k-')
    hold on 
    plot(lWidth.*[1 length(cohTest)],cosieStructCOH.redEML(2,EMLidx).*[1 1],'k-')
    plot(lWidth.*(1:length(cohTest)),cohTest,'r-')
    xlabel('Lateral position (m)')
    ylabel('Coherence Value')
    

   
    subplot(2,2,2)
    plot(lWidth.*[1 length(snrTest)],cosieStructSNR.redEML(1,EMLidx).*[1 1],'k-')
    hold on 
    plot(lWidth.*[1 length(snrTest)],cosieStructSNR.redEML(2,EMLidx).*[1 1],'k-')
    plot(lWidth.*(1:length(snrTest)),snrTest,'b-')
    xlabel('Lateral position (m)')
    ylabel('SNR Value')
    
 
    subplot(2,2,3)
    plot(cosieStructCOH.redEML(1,:),cosieStructCOH.redEML(2,:),'k-.')
    hold on 
    contour(cosieStructCOH.thVector,cosieStructCOH.thVector,abs(1-cosieStructCOH.segSurface),2.^(0.5:-0.5:-10));    
    %xLimsC1 = [min(c1(1,:)) max(c1(1,:))];
    ax = gca;
    xLimsC1 = ax.XLim;
    plot(xLimsC1,xLimsC1,'k-')
    xlim(xLimsC1)
    ylim padded
    xlabel('Lower Threshold')
    ylabel('Upper Threshold')
    axis equal


    subplot(2,2,4)
    plot(cosieStructSNR.redEML(1,:),cosieStructSNR.redEML(2,:),'k-.')
    hold on 
    contour(cosieStructSNR.thVector,cosieStructSNR.thVector,abs(1-cosieStructSNR.segSurface),2.^(0.5:-0.5:-10));
    ax = gca;
    xLimsC2 = ax.XLim;plot(xLimsC2,xLimsC2,'k-')
    xlim(xLimsC2)
    ylim padded
    xlabel('Lower Threshold')
    ylabel('Upper Threshold')
    axis equal

end
