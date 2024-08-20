

function [] = MultipleParametricCohImage(xVals,yVals,rayIdxs,bfImgData,depthVals,powAll,segAll,kLength_BSC_samples,kLength_COH_samples)
    
    %MultipleParametricCohImage(xVals,yVals,rayIdxs,bfImgData,depthVals,powAll,segAll,kLength_BSC_samples,kLength_COH_samples)
    
    B80 = bmode(bfImgData.iq',80);
    
    transpData = zeros(size(bfImgData.fullIM));
    colorData = nan.*zeros(size(bfImgData.fullIM));
    
    xVals = xVals - mean(xVals);

    figure
    ax1 = axes;
    imagesc(ax1,xVals,yVals,B80);
    hold on 
    %plot(ax1,[xVals(rayIdxs(1)) xVals(rayIdxs(end))],yVals(depthIdx).*[1 1],'-','color','red')
    %plot(ax1,[xVals(rayIdxs(end)+5) xVals(rayIdxs(end)+5)], [yVals(axIdxsBSC(1)) yVals(axIdxsBSC(end))],'r^-','MarkerFaceColor','red','LineWidth',2)
    %plot(ax1,[xVals(rayIdxs(1)-5)   xVals(rayIdxs(1)-5)],   [yVals(axIdxsCOH(1)) yVals(axIdxsCOH(end))],'ro-','MarkerFaceColor','red','LineWidth',2)
    
    xlabel('Lateral Position (mm)')
    axis equal
    ylabel('Axial Position (mm)')
    ax1.FontSize = 40;
    ax1.XTick = -16:8:16;
    ax1.YTick = 10:10:50;

    for iDepth = 1:length(depthVals)

        powf0 = powAll(iDepth,:)';
        segBool = segAll(iDepth,:)';

        [~, depthIdx] = min(abs(yVals - depthVals(iDepth)));
        axIdxs = depthIdx-round(kLength_BSC_samples/2) : depthIdx + round(kLength_BSC_samples/2) -1 ;

      
        powfLong = repmat(powf0,1,kLength_BSC_samples);
        hideVectorLong = repmat(segBool,1,kLength_BSC_samples);
        colorData(rayIdxs,axIdxs) = log10(powfLong); 
        transpData(rayIdxs,axIdxs) = 0.6.*hideVectorLong; 

    end
    
    axNew = axes;
    imagesc(axNew,xVals , yVals, colorData','AlphaData',transpData')    
    cB = colorbar;
    cB.Label.String = '\mu (dBcm^{-1}sr^{-1})';
    cB.Label.FontSize= 40;
    cB.FontSize = 30;
    cB.Ticks = round(min(colorData(:))):round(max(colorData(:)));


    linkaxes([ax1,axNew])
    
    axNew.Visible = 'off';
    axNew.XTick = [];
    axNew.YTick = [];
            
  
    %%Give each one its own colormap
    colormap(ax1,'gray')
    axis tight 
    axis equal

    input('Resize')

    axisPosition = get(ax1,'Position');
    colorbarPosition = get(cB, 'Position');

    % Adjust the colorbar height and position
    colorbarPosition(2) = axisPosition(2);  % Match the bottom position
    colorbarPosition(4) = axisPosition(4);  % Match the height
    
    % Apply the new position to the colorbar
    set(cB, 'Position', colorbarPosition);


    cB.Location = 'eastoutside';
    cB.Location = 'manual';

end