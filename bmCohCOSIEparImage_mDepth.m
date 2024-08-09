function [] = bmCohCOSIEparImage_mDepth(xVals,yVals,bfImgData,depthIdx,axIdxsBSC,axIdxsCOH,powf0,segBool,rayIdxs,segEML)

    B80 = bmode(bfImgData.iq',80);

    transpData = zeros(size(bfImgData.fullIM));
    colorData = nan.*zeros(size(bfImgData.fullIM));
    


    figure
    ax1 = axes;
    imagesc(ax1,xVals,yVals,B80);
    hold on 
    plot(ax1,[xVals(rayIdxs(1)) xVals(rayIdxs(end))],yVals(depthIdx).*[1 1],'-','color','red')
    plot(ax1,[xVals(rayIdxs(end)+5) xVals(rayIdxs(end)+5)], [yVals(axIdxsBSC(1)) yVals(axIdxsBSC(end))],'r^-','MarkerFaceColor','red','LineWidth',2)
    plot(ax1,[xVals(rayIdxs(1)-5) xVals(rayIdxs(1)-5)], [yVals(axIdxsCOH(1)) yVals(axIdxsCOH(end))],'ro-','MarkerFaceColor','red','LineWidth',2)
    
    xlabel('Lateral Position (m)')
    axis equal
    ylabel('Axial Position (m)')
     
    
    ax2 = axes;
    powfLong = repmat(powf0,1,size(axIdxsBSC,2));
    hideVectorLong = repmat(segBool,1,size(axIdxsBSC,2));
    colorData(rayIdxs,axIdxsBSC) = log(powfLong); 
    transpData(rayIdxs,axIdxsBSC) = 0.6.*hideVectorLong; 

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

    speckleScore = 100 - (outSeg-segEML);

    title(ax1,['Speckle Score = ',num2str(round(speckleScore)),'%'])
    axis tight 
    axis equal
   

end 

