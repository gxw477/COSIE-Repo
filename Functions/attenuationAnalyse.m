function [aeStruct ] = attenuationAnalyse(bfImgData,slDistCM,endDepthCM,segBool,IDF,rayIdxs,allDepthsCM)
    
    if nargin == 0
        load('C:\Users\gwest\Documents\MATLAB\Temp\data.mat')
    end

    %%%% Input %%%%%%%%%
    % bfImgData : structure produced by beamforming script
    % slDist    : skin to liver distance (cm)
    % endDepth  : maximum depth 
    % segBool   : boolean matrix from segmentation
    % IDF       : inverse diffraction filter from BAE_GUI_GEORGE
    % rayIdxs   : which image lines to use
    % allDepths : depths corresponding to the segBool array
    
    %%%% Output %%%%%%%%
    % aeStructaeStruct = struct;
    % aeStruct.a0    : zero frequency crossing of linear attenuation fit
    % aeStruct.alpha : linear frequency dependence of linear attenuation
    % aeStruct.corr  : RMSE error 
    % aeStruct.fL    : lower frequency of fit
    % aeStruct.fU    : upper frequency of fit
    
    % aeStruct.a0_filt    : a0 but with segmentation data
    % aeStruct.alpha_filt : alpha but for segmentation data
    % aeStruct.corr_filt  : corr but for segmentation data



    %Window size
    wSize = IDF.Filt.samp_size;
    %number of frequencies
    nF = wSize;
    %conversion sample 2 cm 
    samp2cm= bfImgData.lambda/(4*2)*1e2;
    %window size in cm
    wSizeCM = wSize*samp2cm/2;
    %Window positions in cm
    IDF.Filt.wPosCM = IDF.Filt.wpos*samp2cm;
    
    %frequencies
    f = (0:(wSize-1)).*(bfImgData.fs/(wSize));
    
    %IDFfilt = IDF.Filt.F;
    
    %Hann window function
    win = sqrt(8/3).*hann(wSize);

    arrayDepthsLiver = [(allDepthsCM(1)-wSizeCM/2),slDistCM + wSizeCM/2, endDepthCM+wSizeCM/2];
    arrayDepthsIDF = [(IDF.Filt.wpos(1).*samp2cm - wSizeCM/2), (IDF.Filt.wpos(end).*samp2cm +wSizeCM/2)];

    %fast time indices within the range of analysis and inside the liver 
    idxBoolLiver = bfImgData.yVals.*1e2 > min(arrayDepthsLiver) & bfImgData.yVals.*1e2 < max(arrayDepthsLiver);
    
    %fast time indices in the range of the IDF filter 
    idxBoolIDF = bfImgData.yVals.*1e2 > arrayDepthsIDF(1) &  bfImgData.yVals.*1e2 < arrayDepthsIDF(2); 

    idxBool = idxBoolLiver & idxBoolIDF;


    %Start fast time index
    startIdxRF = find(idxBool,1);
    endIdxRF   = find(idxBool,1,'last');
    %fast time posns in cm
    zVals = bfImgData.yVals(startIdxRF:endIdxRF).*1e2;

    %window position indices 
    startIdxWin = startIdxRF + wSize/2;
    endIdxWin = endIdxRF -wSize/2;

    %centre position of windows used
    winZvalsCM = bfImgData.yVals(startIdxWin:endIdxWin).*1e2;
   

    %number of window positions to work with
    nDepthWin = endIdxWin - startIdxWin + 1;

    %number of raylines 
    nRays = length(rayIdxs);    
        
    %interpolate segBool to match the window spacing for spect analysis
    %interpolate segBool to match RF data dimensions
    [xQ1,zQ1] = meshgrid(1:nRays,allDepthsCM);
    [xQ2,zQ2] = meshgrid(1:nRays,winZvalsCM);
    segBool2 = interp2(xQ1,zQ1, double(segBool'), xQ2, zQ2 ); 
    segBool2(isnan(segBool2)) = 0;
    segBool2 = round(segBool2');

    figure
    subplot(1,2,1)
    imagesc(1:96,allDepthsCM,segBool')
    subplot(1,2,2)
    imagesc(1:96,winZvalsCM,round(segBool2'))
    close all

    %power spectrum variables for the backscattered data
    spect = zeros(nRays,nDepthWin,wSize/2);
    spectCorr = spect;
    spectCorrFilt = spectCorr;
    
    %zvals for fitting in cm
    zAtt = repmat(winZvalsCM,[length(rayIdxs),1]);
    zAtt_mask = nan(size(zAtt)); 

    %Query values from IDF 
    [zQ,fQ1] = meshgrid(IDF.Filt.wPosCM,IDF.Filt.f);
    %Query values from yVals 
    [yQ,fQ2] = meshgrid(bfImgData.yVals.*1e2,IDF.Filt.f);

    %interpolate IDF to match yVals
    
    idfInterp = interp2(zQ,fQ1,IDF.Filt.F,yQ,fQ2);
 
    
    for iDepth = 1:nDepthWin
        
        iDepth;

        if iDepth == 132
            iDepth;
        end

        %fast time indices for sampling
        depthIdxs = startIdxRF + (iDepth-1) + (0:wSize-1);

        idfCorrSpect = idfInterp(:,round(mean(depthIdxs)))';
              
        for iRay = 1:nRays
            
            %calculate spectrum
            temp = abs(fft(bfImgData.fullIM(rayIdxs(iRay),depthIdxs).* win')).^2;
            %raw
            spect(iRay,iDepth,:) = temp(1:wSize/2);
            %IDF corrected 
            spectCorr(iRay,iDepth,:) = temp(1:wSize/2)./idfCorrSpect;

            %fill if segBool2 is true 
            if segBool2(iRay,iDepth)
                %Unsegmented value
                spectCorrFilt(iRay,iDepth,:) = spectCorr(iRay,iDepth,:);
                zAtt_mask(iRay,iDepth) = zAtt(iRay,iDepth);
            end
        end
    end

    [~,fL] = min(abs(4.5e6- f));
    [~,fU] = min(abs(5.7e6 - f));

    a = zeros(1,nF);
    b = a;
    a2 = b;
    b2 = a2;

    for iF = fL:fU      
        
        %results are dB , dB/cm, [], y(iF) = a(iF)+b(iF)*2*x
        
        %log of data
        Flog = mean(10*log10(spectCorr(:,:,iF)),1)';

        zeroBool = spectCorrFilt(:,:,iF) == 0;

        for iRayBool = 1:size(zeroBool,1)
            spectCorrFilt(iRayBool,zeroBool(iRayBool,:),iF) = nan;
        end

        FlogFilt = mean(10*log10(spectCorrFilt(:,:,iF)),1,'omitnan');
    
        options1 = [2*0.5*(endDepthCM + slDistCM) 0 0];
        
        zFit1 = 2*zAtt;

        %zAtt_mask should be either nan or the identical values, so we can
        %do a mean to reduce dimensionality and remove the nans
        zFit2 = 2*mean(zAtt_mask,1,'omitnan');
        yFit1 = Flog;
        yFit2 = FlogFilt;
        
        %results in dB/cm
        %[a(iF),b(iF),c(iF),e(iF,:),re(iF)] = lsqn(zFit1,yFit1,options1);
        coeffs = polyfit(zFit1(1,:),yFit1(:),1);
        a(iF) = coeffs(2);
        b(iF) = coeffs(1);
    
        validIdx  = ~isnan(zFit2) & ~isnan(yFit2);
        zFit2 = zFit2(validIdx);
        
        %bool5 = zFit2 > 5 & zFit2 < 5.008;

        yFit2 = yFit2(validIdx);
    
        coeffs = polyfit(zFit2(validIdx),yFit2(validIdx),1);
        a2(iF) = coeffs(2);
        b2(iF) = coeffs(1);
    
        %[a2(iF),b2(iF),c2(iF),e2(iF,:),re2(iF)] = lsqn(zFit2,yFit2,options1);
        
    
        if iF==26
            figure 
            subplot(1,4,1)
            plot(zFit1(1,:),yFit1(:),'k.','MarkerFaceColor','k')
            hold on 
            plot(zFit1(1,:),a(iF)+b(iF).*zFit1(1,:),'k-.','LineWidth',2)
            plot(zFit2,yFit2(:),'r.')
            plot(zFit2,a2(iF)+b2(iF).*zFit2,'r-*','LineWidth',2)
            
        end
    
    end


    b  = -b; 
    b2 = -b2;

    [~,iF0 ] = min(abs(f - bfImgData.Trans.frequency*1e6));

    %results are dB/cm, dB/cm/MHz
    temp  = polyfit(f(fL:fU)./1e6,b(fL:fU),1);
    temp2 = polyfit(f(fL:fU)./1e6,b2(fL:fU),1);
    
    %a0 = temp(2);
    %alpha = temp(1);
    %a02 =   temp2(2);
    %alpha2 = temp2(1);
    
    [a0,  alpha  ] = lsqn(f(fL:fU)./1e6, b(fL:fU) ,[f(iF0)./1e6 (f(2)-f(1))./1e6 2]); 
    [a02, alpha2 ] = lsqn(f(fL:fU)./1e6, b2(fL:fU),[f(iF0)./1e6 (f(2)-f(1))./1e6 2]); 
    
    if 1
        subplot(1,4,2)
        plot(f(fL:fU),b(fL:fU),'ko','MarkerFaceColor','k')
        hold on
        plot(f(fL:fU), a0+ alpha.*f(fL:fU)./1e6,'k-.')
        plot(f(fL:fU), b2(fL:fU),'ro','MarkerFaceColor','r')
        plot(f(fL:fU), a02 + alpha2*f(fL:fU)./1e6,'r-.')
        title(['End Depth : ',num2str(endDepthCM)])
      
        subplot(1,4,3)
        imagesc(1:96,2.*winZvalsCM,segBool2')
        
     
        subplot(1,4,4)
        imagesc(1:96,2.*allDepthsCM,segBool')
        ylim(2.*[min(winZvalsCM) max(winZvalsCM)])
        close all


    end

    resid =      a0  + alpha.*f(fL:fU)./1e6 - b(fL:fU);
    resid_Filt = a02 + alpha2.*f(fL:fU)./1e6 - b2(fL:fU);

    aeStruct = struct;
    aeStruct.a0 = a0;
    aeStruct.alpha = alpha; 
    aeStruct.corr = sqrt(mean(resid.^2 , 2,'omitnan'));
    aeStruct.fL = fL;
    aeStruct.fU = fU;
    
    aeStruct.a0_filt = a02;
    aeStruct.alpha_filt = alpha2; 
    aeStruct.corr_filt = sqrt(mean(resid_Filt.^2 , 2,'omitnan'));
    aeStruct.fL = fL;
    aeStruct.fU = fU;
    

end