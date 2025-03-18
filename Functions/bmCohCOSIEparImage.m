function [] = bmCohCOSIEparImage(xVals,yVals,bfImgData,depthIdx,axIdxsBSC,axIdxsCOH,powf0,segBool,rayIdxs,segEML,titleString)
% bmCohCOSIEparImage(xVals,yVals,bfImgData,depthIdx,axIdxsBSC,axIdxsCOH,powf0,segBool,rayIdxs,segEML,titleString)
%
%   I 
%   ~~~~~~~~~~
%   xVals & yVals vector (1xM lines, 1:k samples) (mm)!!! 
%   bfImgData     struct
%   depthIdx    scalar 
%   axIdxsBSC   vector (1xN idxs) axial indices for BSC estm.
%   axIdxsCOH   vector (1xN idxs)    "      "  for Coh estm. 
%   powf0       vector (1xM lines) absolute backscattered power at f0
%   segBool     keep = 1, seg = 0
%   rayIdxs     vector (1xM' lines) lines to keep
%   segEML      matrix (2xS seg points) 
    
    %rawIQ{1} = bfImgData.IData{1}(:,:,1) +1i.*bfImgData.QData{1}(:,:,1);
    rawIQ{1} = bfImgData.iq';

    [iqV , ~ , zOut] = demodulateIQfn(bfImgData.PData,rawIQ,2);
    nSamps = length(zOut);
    bModeViq = bmode(iqV,70);

    transpData = zeros(size(bfImgData.fullIM));
    colorData = nan.*zeros(size(bfImgData.fullIM));
    
    xVals = xVals - mean(xVals);
    

    bmXVals = 1e3.*(bfImgData.PData.Origin(1) + (1:size(bModeViq,2)).*bfImgData.PData.PDelta(1))*bfImgData.lambda;
    %bmYVals = (bfImgData.PData.Origin(3) + (1:size(bModeViq,1)).*bfImgData.PData.PDelta(3))*bfImgData.lambda;
   


    %bmZvals = 1e3.*bfImgData.lambda.*(bfImgData.PData.Origin(3)+ zOut(1:nSamps).*bfImgData.PData.PDelta(3));
    
    
    figure
    ax1 = axes;
    imagesc(ax1,xVals,yVals,bModeViq);
    hold on 
    plot(ax1,[xVals(rayIdxs(1)) xVals(rayIdxs(end))],yVals(depthIdx).*[1 1],'-','color','red')
    plot(ax1,[xVals(rayIdxs(end)+5) xVals(rayIdxs(end)+5)], [yVals(axIdxsBSC(1)) yVals(axIdxsBSC(end))],'r^-','MarkerFaceColor','red','LineWidth',2)
    plot(ax1,[xVals(rayIdxs(1)-5) xVals(rayIdxs(1)-5)], [yVals(axIdxsCOH(1)) yVals(axIdxsCOH(end))],'ro-','MarkerFaceColor','red','LineWidth',2)
    
    xlabel('Lateral Position (mm)')
    %axis equal
    ylabel('Axial Position (mm)')
     
    
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

    %title(ax1,[titleString,' Speckle Score = ',num2str(round(speckleScore))])
    axis tight 
    %axis equal

    %input('Resize : ')

    axisPosition = get(ax1,'Position');
    colorbarPosition = get(cB, 'Position');

    % Adjust the colorbar height and position
    colorbarPosition(2) = axisPosition(2);  % Match the bottom position
    colorbarPosition(4) = axisPosition(4);  % Match the height
    
    % Apply the new position to the colorbar
    set(cB, 'Position', colorbarPosition);


    cB.Location = 'eastoutside';
    cB.Location = 'manual';

    %ylim([5 45])
    %xlim([-17 17])
    axis tight

end 

