

function [] = MultipleParametricCohImage(xVals,yVals,rayIdxs,bfImgData,depthVals,powAll,segAll,kLength_BSC_samples,kLength_COH_samples)
    
    %MultipleParametricCohImage(xVals,yVals,rayIdxs,bfImgData,depthVals,powAll,segAll,kLength_BSC_samples,kLength_COH_samples)
    
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
    
    
    for iDepth = 1
        
        axNew(i) = axes;
        powfLong = repmat(powf0,1,kLength_BSC_samples);
        hideVectorLong = repmat(segBool,1,kLength_BSC_samples);
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
            
    end
    
    
  
    %%Give each one its own colormap
    colormap(ax1,'gray')


end