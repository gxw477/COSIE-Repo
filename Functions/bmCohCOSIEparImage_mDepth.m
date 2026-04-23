function [ax1,ax2,cB] = bmCohCOSIEparImage_mDepth(xVals,yVals,bfImgData,iqBool,axIdxsBSC,kLength,powf0,segBool,rayIdxs,segEML,titleString)
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
    
    if iqBool
        
        Itemp = bfImgData.IData{1}(:,:,1);
        Qtemp = bfImgData.QData{1}(:,:,1);
        PData = bfImgData.PData;

        %axial pix sep'n 
        dz = PData.PDelta(3); 

        %number of lines
        nr = size(Itemp,2);
        %number of time samples
        l  = size(Itemp,1);

        interpFactor = 8;
        phasor=repmat([0:l-1]'*2*pi*dz,1,nr);
        phasor=reshape(phasor,l,nr);
        %IQ=interp1([0:l-1],IQData{1}(1:l,:,:),[0:0.25/dz:l-1],'cubic');
        IQa=interp1([0:l-1],abs(Itemp),[0:(dz/interpFactor):l-1],'cubic');
        IQp=interp1([0:l-1],unwrap(angle(Qtemp)-phasor),[0:(dz/interpFactor):l-1],'cubic');
        IQ=IQa.*exp(sqrt(-1)*IQp);        
        bModeViq = bmode(IQ,80);
        
        bmXvals = ((0:nr-1)).*bfImgData.lambda*PData.PDelta(1)*1e3;
        bmZvals = (PData.Origin(3)+ (0:dz/interpFactor:l-1)).*bfImgData.lambda*PData.PDelta(3)*1e3;

    else

        bModeViq = bmode(bfImgData.iq',75);

        bmXvals = (1:size(bfImgData.fullIM,1)).*bfImgData.dX.*1e3;
        bmZvals = 1e3.*bfImgData.lambda.*(bfImgData.PData.Origin(3)+ (1:2*bfImgData.PData.Size-1).*0.5*bfImgData.PData.PDelta(3));

    end

    bmXvals = bmXvals - mean(bmXvals);
        
    colorData = nan.*zeros(size(xVals,2),size(bfImgData.fullIM,2));
    transpData = zeros(size(colorData));

    xVals = xVals - mean(xVals);
        
    %bModeViq2 = sigmoidEnhance(bModeViq,5,0.5);
    bModeViq2 = interp2(double(bModeViq),2);%claheEnhance(bModeViq,'Distribution', 'rayleigh', 'Alpha', 0.1,'ClipLimit',0.5);
  
    figure
    imagesc(bmXvals,bmZvals,bModeViq2);
    xlabel('Lateral Position (mm)')
    ylabel('Axial Position (mm)')
    axis equal
    axis tight 
    ylim([5 55])

    figure
    ax1 = axes;
    imagesc(ax1,bmXvals,bmZvals,bModeViq2);   
    xlabel('Lateral Position (mm)')
    ylabel('Axial Position (mm)')
    axis equal
    
    ax2 = axes;
    powfLong = repelem(powf0,1,kLength);
    hideVectorLong = repelem(segBool,1,kLength);
    colorData(rayIdxs,axIdxsBSC) = log10(powfLong); 
    transpData(rayIdxs,axIdxsBSC) = 0.5.*hideVectorLong; 

    imagesc(ax2,xVals , yVals, colorData','AlphaData',transpData')    
    cB = colorbar;
    cB.Label.String = 'log_{10}(\kappa)';
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

