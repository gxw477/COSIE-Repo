function [] = depthFeaturePlot(cosieStructCOH,cosieStructSNR,EMLidx,lWidth,cohTest,snrTest,newfigBool)

    % depthFeaturePlot(cosieStructCOH,cosieStructSNR,EMLidx,lWidth,cohTest,snrTest)
    % cosieStructCOH : cosie struct containing EML for coherence seg
    % cosieStructSNR : cosie struct containing EML for snr seg
    % EMLidx : self explanatory
    % lWidth : image line width
    % cohTest : coherence values in image 
    % snrTest : snr values for image


    
    figure 
    subplot(1,2,1)
    plot(lWidth.*[1 length(cohTest)],cosieStructCOH.redEML(1,EMLidx).*[1 1],'k-')
    hold on 
    plot(lWidth.*[1 length(cohTest)],cosieStructCOH.redEML(2,EMLidx).*[1 1],'k-')
    plot(lWidth.*(1:length(cohTest)),cohTest,'r-')
    xlabel('Lateral position (m)')
    ylabel('Coherence Value')
    

   
    subplot(1,2,2)
    plot(lWidth.*[1 length(snrTest)],cosieStructSNR.redEML(1,EMLidx).*[1 1],'k-')
    hold on 
    plot(lWidth.*[1 length(snrTest)],cosieStructSNR.redEML(2,EMLidx).*[1 1],'k-')
    plot(lWidth.*(1:length(snrTest)),snrTest,'b-')
    xlabel('Lateral position (m)')
    ylabel('SNR Value')
    
    

end
