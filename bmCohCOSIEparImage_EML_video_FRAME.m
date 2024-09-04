function [] = bmCohCOSIEparImage_EML_video_FRAME(xVals,yVals,bfImgData,depthIdx,axIdxs,powf0,segBool,rayIdxs,speckleCOSIE,cohTest,sumIdx,segIdx)

    %EML is the value of the EML up until the current index (for plotting
    %purposes only)
    

    B80 = bmode(bfImgData.iq',80);

    transpData = zeros(size(bfImgData.fullIM));
    colorData = nan.*zeros(size(bfImgData.fullIM));
    
    
    

    figure
    subplot(1,2,1)
    ax1 = gca;
    imagesc(ax1,xVals,yVals,B80);
    hold on 
    plot(ax1,[xVals(rayIdxs(1)) xVals(rayIdxs(end))],yVals(depthIdx).*[1 1],'-','color','red')
    plot(ax1,[xVals(rayIdxs(end)) xVals(rayIdxs(end))], [yVals(axIdxs(1)) yVals(axIdxs(end))],'-','color','red')
    xlabel('Lateral Position (m)')
    axis equal
    ylabel('Axial Position (m)')
     
    
    ax2 = axes;
    ax2.Position = ax1.Position 
    powfLong = repmat(powf0,1,size(axIdxs,2));
    
    hideVectorLong = repmat(segBool,1,size(axIdxs,2));
    colorData(rayIdxs,axIdxs) = log(powfLong); 
    transpData(rayIdxs,axIdxs) = 0.6.*hideVectorLong; 

    imagesc(ax2,xVals , yVals, colorData','AlphaData',transpData')    
    cB = colorbar;
    cB.Label.String = 'Log(\mu)';
    cB.Label.FontSize= 10;
    cB.FontSize = 8;
    cB.Position = [0.4432 0.1835 0.0300 0.6663];
    %cB.Limits = [12 20];
    
    linkaxes([ax1,ax2])
    
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    %%Give each one its own colormap
    colormap(ax1,'gray')
    %colormap(ax2,'winter')
    %ylim([0.04 0.09])
    
    outSeg = 100*(1-length(find(segBool))/length(segBool));

    title(ax1,['Seg = ',num2str(outSeg),'%'])
    axis tight 
    axis equal
   
    [speckleCOHpdf,speckleCOHbins] = histcounts(speckleCOSIE.cohSum);

    speckleCOHbins = speckleCOHbins(1:end-1)+(speckleCOHbins(2)-speckleCOHbins(1))/2;
    speckleCOHpdf = speckleCOHpdf./(sum(speckleCOHpdf));
    
     
    
    subplot(2,2,4)
   
    hold on  
    set(gca,'YColor','k')
    yticks([])
    yyaxis right 
    set(gca,'YColor','k')
    h = histogram(cohTest,'Normalization','pdf','BinWidth',2);
    globMax = max([speckleCOHpdf,h.Values]);
    plot([speckleCOSIE.EML(1,segIdx),speckleCOSIE.EML(1,segIdx)],[0 1.2*globMax],'r-','LineWidth',0.5)
    plot([speckleCOSIE.EML(2,segIdx),speckleCOSIE.EML(2,segIdx)],[0 1.2*globMax],'r-','LineWidth',0.5)
    plot(speckleCOHbins,speckleCOHpdf,'k-','LineWidth',2)
    xlabel('Coherence Value')
    ylabel('PDF')
    ylim([0 1.2*globMax])

    subplot(2,2,2)
    set(gca,'YColor','k')
    yticks([])
    yyaxis right 
    set(gca,'YColor','k')
    plot(speckleCOSIE.EML(1,1:segIdx),speckleCOSIE.EML(2,1:segIdx),'k-o','MarkerFaceColor','k','MarkerSize',2)
    hold on 
    plot([-sumIdx sumIdx],[-sumIdx sumIdx],'k-')
    contour(speckleCOSIE.thVector,speckleCOSIE.thVector,abs(1-speckleCOSIE.segSurface),2.^(0.5:-0.5:-10))
    cB2 = colorbar;
    cB2.Location = 'westtoutside';
    cB2.Label.String = 'Seg. Strength';
    cB2.Label.FontSize = 9;
    xlim([-2 max(speckleCOSIE.thVector)])
    ylim([-2 max(speckleCOSIE.thVector)])
    xlabel('Lower Threshold')
    ylabel('Upper Threshold')
    %axis equal

end 

