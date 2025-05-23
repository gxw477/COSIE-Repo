function [ax1,ax2,cB] = bmCohCOSIEparImage_mDepth_C16D(xVals,yVals,bfImgData,depthIdx,axIdxsBSC,kLength,powf0,segBool,rayIdxs,segEML,titleString)
% bmCohCOSIEparImage_mDepth(xVals,yVals,bfImgData,depthIdx,axIdxsBSC,axIdxsCOH,powf0,segBool,rayIdxs,segEML,titleString)
%
%   I 
%   ~~~~~~~~~~
%   xVals & yVals vector (1xM lines, 1:k samples) (mm)!!! 
%   bfImgData     struct
%   depthIdx    scalar (only plots the hzal line now 
%   axIdxsBSC   vector (DxN idxs) axial indices for BSC estm.
%   kLength     kernel length
%   powf0       vector (DxM lines) absolute backscattered power at f0
%   segBool     keep = 1, seg = 0
%   rayIdxs     vector (1xM' lines) lines to keep
%   segEML      matrix (2xS seg points) 
    
    
    iqV = rf2iq(bfImgData.fullIM,bfImgData.fs);
    bModeViq = bmode(iqV',100);

    transpData = zeros(size(bfImgData.fullIM));
    colorData = nan.*zeros(size(bfImgData.fullIM));
    
    xVals = xVals - mean(xVals);
    

    bmXVals = 1e3.*(bfImgData.PData.Origin(1) + (1:size(bModeViq,2)).*bfImgData.PData.PDelta(1))*bfImgData.lambda;
    bmZvals = 1e3.*bfImgData.lambda.*(bfImgData.PData.Origin(3)+ (1:2*bfImgData.PData.Size-1).*0.5*bfImgData.PData.PDelta(3));

    

    figure
    ax1 = axes;
    imagesc(ax1,xVals,yVals,interp2(double(bModeViq),2));
    hold on 
    %plot(ax1,[xVals(rayIdxs(1)) xVals(rayIdxs(end))],yVals(depthIdx).*[1 1],'-','color','red')
    %plot(ax1,[xVals(rayIdxs(end)+5) xVals(rayIdxs(end)+5)], 10.*[1 1],'r^-','MarkerFaceColor','red','LineWidth',2)
    %plot(ax1,[xVals(rayIdxs(1)-5) xVals(rayIdxs(1)-5)], [yVals(axIdxsCOH(1)) yVals(axIdxsCOH(end))],'ro-','MarkerFaceColor','red','LineWidth',2)
    
    xlabel('Lateral Position (mm)')
    axis equal
    ylabel('Axial Position (mm)')
     
    
    ax2 = axes;
    powfLong = repelem(powf0,1,kLength);
    hideVectorLong = repelem(segBool,1,kLength);
    colorData(rayIdxs,axIdxsBSC) = log(powfLong); 
    transpData(rayIdxs,axIdxsBSC) = 0.5.*hideVectorLong; 

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

    %title(ax1,['Speckle Score = ',num2str(round(speckleScore)),'%'])
    axis tight 
    axis equal
   
    title(titleString)
end 

