function [] = EMLPlotter(speckleCOSIE)

    figure
    subplot(1,2,1)
    plot(speckleCOSIE.EML(1,:),speckleCOSIE.EML(2,:),'k-.')
    hold on 
    contour(speckleCOSIE.thVector,speckleCOSIE.thVector,abs(1-speckleCOSIE.segSurface),2.^(0.5:-0.5:-10))
    plot([-10 20],[-10 20],'k-')
    xlim([min(speckleCOSIE.EML(1,:))-1 max(speckleCOSIE.EML(1,:))+1])
    ylim padded
    
    xlabel('Lower Threshold')
    ylabel('Upper Threshold')
    axis equal

    subplot(2,2,2)
    histogram(speckleCOSIE.cohSum,'Normalization','pdf')
    xlabel('Coherence Value')
    ylabel('P.D.F.')

    subplot(2,2,4)
    plot(speckleCOSIE.cohSum,speckleCOSIE.powf0,'k.')
    xlabel('Coherence')
    ylabel('Power (a.u.)')
    cc = corrcoef(speckleCOSIE.cohSum,speckleCOSIE.powf0);
    title(['R = ',num2str(cc(2))])

end