function [] = bmCohCOSIEparImage(xVals,yVals,bfImgData,depthIdx,axIdxs,powf0,segBool,rayIdxs)

    B80 = bmode(bfImgData.iq',80);

    transpData = zeros(size(bfImgData.fullIM));
    colorData = nan.*zeros(size(bfImgData.fullIM));
    
    
    ax1 = axes;
    imagesc(ax1,xVals,yVals,B80);
    hold on 
    plot(ax1,[xVals(rayIdxs(1)) xVals(rayIdxs(end))],yVals(depthIdx).*[1 1],'-','color','red')
    plot(ax1,[xVals(rayIdxs(end)) xVals(rayIdxs(end))], [yVals(axIdxs(1)) yVals(axIdxs(end))],'-','color','red')
    xlabel('Lateral Position (m)')
    axis equal
    ylabel('Axial Position (m)')
     
    
    ax2 = axes;
    powfLong = repmat(powf0,1,size(axIdxs,2));
    
    hideVectorLong = repmat(segBool,1,size(axIdxs,2));
    colorData(rayIdxs,axIdxs) = log(powfLong); 
    transpData(rayIdxs,axIdxs) = 0.6.*hideVectorLong; 

    imagesc(ax2,xVals , yVals, colorData','AlphaData',transpData')    
    cB = colorbar;
    cB.Label.String = 'Log(\kappa)';
    cB.Label.FontSize= 20;
    cB.Position(3) = 0.02;
    
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
   


end 

