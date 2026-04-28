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
    %aeStruct.a0    : zero frequency crossing of linear attenuation fit
    %aeStruct.alpha : linear frequency dependence of linear attenuation
    %aeStruct.corr  : RMSE error 
    %aeStruct.fL    : lower frequency of fit
    %aeStruct.fU    : upper frequency of fit
    
    %aeStruct.a0_filt    : a0 but with segmentation data
    %aeStruct.alpha_filt : alpha but for segmentation data
    %aeStruct.corr_filt  : corr but for segmentation data



    %Window size
    wSize = 2*size(IDF.Filt.F,1);
    %number of frequencies
    nF = wSize;
    %conversion sample 2 cm 
    samp2cm= bfImgData.lambda/(4*2)*1e2;
    %window size in cm
    wSizeCM = wSize*samp2cm/2;
    %Window positions in cm
    wPosCM = IDF.Filt.wpos*samp2cm;
    
    %frequencies
    f = (0:(wSize-1)).*(bfImgData.fs/(wSize));
    
    %IDFfilt = IDF.Filt.F;
    
    %Hann window function
    win = sqrt(8/3).*hann(wSize);

    %fast time indices within the range of analysis and inside the liver 
    idxBoolLiver = bfImgData.yVals.*1e2 >= (allDepthsCM(1)-wSizeCM/2) & bfImgData.yVals.*1e2 >= slDistCM + wSizeCM/2 & bfImgData.yVals.*1e2 <= endDepthCM+wSizeCM/2;
    
    %fast time indices in the range of the IDF filter 
    idxBoolIDF = bfImgData.yVals.*1e2 >= (IDF.Filt.wpos(1).*samp2cm - wSizeCM/2) &  bfImgData.yVals.*1e2 <= (IDF.Filt.wpos(end).*samp2cm +wSizeCM/2); 

    idxBool = idxBoolLiver & idxBoolIDF;

    %zVals available based on idxBool
    zVals = bfImgData.yVals(idxBool);



    %ensure that the lengths are divisible  (might not be due to placement of segmentation windows and fast time sampling)
    while round(length(zVals)/length(allDepthsCM)) ~= length(zVals)/length(allDepthsCM)
        zVals = zVals(1:end-1);
    end

    interpRatio = length(zVals)/length(allDepthsCM);
   
    %[segBoolRep] = interp2(1:size(segBool,1),segBoolZVals./100,  )
    segBool2 = false(length(rayIdxs),length(zVals));

    for i = 1:length(rayIdxs)
        segBool2(i,:)= repelem(segBool(i,:),interpRatio);
    end
    
    %Start fast time index
    startIdxRF = find(idxBool,1);
    endIdxRF   = find(idxBool,1,'last');
  
    %window position indices 
    startIdxWin = startIdxRF + wSize/2;
    endIdxWin = endIdxRF -wSize/2;
    
    %fast time posns in cm
    winZvalsCM = bfImgData.yVals(startIdxWin:endIdxWin).*1e2;

    

    %number of window positions to work with
    nDepthWin = endIdxWin - startIdxWin + 1;
    %number of fast time axis 
    nDepthTime = length(zVals);

    %number of raylines 
    nRays = length(rayIdxs);    

    %power spectrum variables for the backscattered data
    spect = zeros(nRays,nDepthWin,wSize/2);
    spectCorr = spect;
    spectCorrFilt = spectCorr;
    
    %zvals for fitting in cm
    zAtt = repmat(winZvalsCM,[length(rayIdxs),1]);
    zAtt_mask = nan(size(zAtt));

    %match IDF to the image sampling

    for iDepth = 1:nDepthWin
        
        iDepth;

        %fast time indices for sampling
        depthIdxs = startIdxRF + (iDepth-1) + (0:wSize-1);

        idfCorrSpect = interp2(IDF.Filt.f,wPosCM,IDF.Filt.F',IDF.Filt.f,zVals(iDepth)*1e2);

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

    [~,fL] = min(abs(4.8e6- f));
    [~,fU] = min(abs(5.7e6 - f));

    a = zeros(1,nF);
    b = a;
    a2 = b;
    b2 = a2;

    for iF = fL:fU      
        
        %results are dB , dB/cm, [], y(iF) = a(iF)+b(iF)*2*x
        
        %log of data
        Flog = 10*log10(spectCorr(:,:,iF));
        FlogFilt = 10*log10(spectCorrFilt(:,:,iF));
    
        %options1 = [2*0.5*(endDepth + slDist) 0 0];
        
        zFit1 = 2*zAtt;
        zFit2 = 2*zAtt_mask;
        yFit1 = Flog;
        yFit2 = FlogFilt;
        
        %results in dB/cm
        coeffs = polyfit(zFit1(:),yFit1(:),1);
        a(iF) = coeffs(2);
        b(iF) = coeffs(1);
    
        validIdx  = ~isnan(zFit2) & ~isnan(yFit2);
    
        coeffs = polyfit(zFit2(validIdx),yFit2(validIdx),1);
        a2(iF) = coeffs(2);
        b2(iF) = coeffs(1);
    
        %[a2(iF),b2(iF),c2(iF),e2(iF,:),re2(iF)] = lsqn(xFit2,yFit2,options1);
        
    
        if iF==26
            figure 
            subplot(1,2,1)
            plot(zFit1(:),yFit1(:),'k.','MarkerFaceColor','k')
            hold on 
            plot(zFit1(:),a(iF)+b(iF).*zFit1(:),'k-.','LineWidth',2)
            plot(zFit2(:),yFit2(:),'r.')
            plot(zFit2(validIdx),a2(iF)+b2(iF).*zFit2(validIdx),'r-*','LineWidth',2)
            
        end
    
    end


    b  = -b; 
    b2 = -b2;

    [~,iF0 ] = min(abs(f - bfImgData.Trans.frequency*1e6));

    %results are dB/cm, dB/cm/MHz
    temp  = polyfit(f(fL:fU)./1e6,b(fL:fU),1);
    temp2 = polyfit(f(fL:fU)./1e6,b2(fL:fU),1);
    
    a0 = temp(2);
    alpha = temp(1);
    a02 =   temp2(2);
    alpha2 = temp2(1);
    
    %[a0,  alpha  ] = lsqn(f(fL:fU)./1e6, b(fL:fU) ,[f(iF0)./1e6 (f(2)-f(1))./1e6 2]); 
    %[a02, alpha2 ] = lsqn(f(fL:fU)./1e6, b2(fL:fU),[f(iF0)./1e6 (f(2)-f(1))./1e6 2]); 
    
    if 1
        subplot(1,2,2)
        plot(f(fL:fU),b(fL:fU),'ko','MarkerFaceColor','k')
        hold on
        plot(f(fL:fU), a0+ alpha.*f(fL:fU)./1e6,'k-.')
        plot(f(fL:fU), b2(fL:fU),'ro','MarkerFaceColor','r')
        plot(f(fL:fU), a02 + alpha2*f(fL:fU)./1e6,'r-.')
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