function [aeStruct ] = attenuationAnalyse(bfImgData,slDist,endDepth,segBool,IDF,rayIdxs,allDepths)

    %attenuationAnalyse(BAE,slDist,endDepth,segBool) 
    % BAE - output from BAE_GUI 
    % slDist - start depth (start of liver) cm
    % endDepth - depth of interest  cm 
    % segBool - segmentation    
    
    if nargin == 0
        load('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\test.mat')
    end

    wSize = 2*size(IDF.Filt.F,1);
    samp2cm= bfImgData.lambda/(4*2)*1e2;
    wPosCM = IDF.Filt.wpos*samp2cm;
    wSizeCM = wSize*samp2cm;

    f =  (0:(wSize-1)).*(bfImgData.fs/(wSize));

    IDFfilt = IDF.Filt.F;
    IDFfilt = IDFfilt./max(IDFfilt(:));
    
    win = sqrt(8/3).*hann(wSize);

    nStartwPos = find(wPosCM>=slDist+ wSizeCM/2,1,'first');
    nEndwPos =   find(wPosCM>=endDepth + wSizeCM/2,1,'first');

    nDepth = nEndwPos - nStartwPos + 1;
        
    wPosCM2 = wPosCM(nStartwPos:nEndwPos);

    spect = zeros(nDepth,length(rayIdxs),wSize/2);
    spectCorr = spect;

    spectCorrFilt = nan(size(spectCorr));

    segBoolAtt = false(length(wPosCM),length(rayIdxs));

    for iDepth2 = 1:length(allDepths)
        
        [~,idxs] = find(abs(wPosCM + wSizeCM/2 - allDepths(iDepth2)/10) < 0.5/2);

        boolTemp = segBool(:,iDepth2);

        segBoolAtt(idxs,:) = repmat(boolTemp,[1,length(idxs)])';

    end

    xAtt = repmat(wPosCM2',[1,length(rayIdxs)]);

    xAtt_mask = nan(size(xAtt));

    for iDepth = 1:nDepth

        depthIdx = iDepth + nStartwPos-1;

        yBool = abs(bfImgData.yVals*1e2-wPosCM2(iDepth))<=wSize*samp2cm/2;
        
        if length(find(yBool)) == 79 
            foo = find(yBool);
            yBool(foo(end)+1) = true;
        end

        for iLine = 1:length(rayIdxs)
            temp = abs(fft(bfImgData.fullIM(rayIdxs(iLine),yBool).* win')).^2;
            spect(iDepth,iLine,:) = temp(1:wSize/2);
            spectCorr(iDepth,iLine,:) = temp(1:wSize/2)./IDFfilt(:,depthIdx)';
        
            
            if (segBoolAtt(iDepth,iLine))
                xAtt_mask(iDepth,iLine) = xAtt(iDepth,iLine);
                spectCorrFilt(iDepth,iLine,:) = squeeze(spectCorr(iDepth,iLine,:));

          
            end
        
        
        end


    end
    


    xAtt2 = xAtt;
    xAtt_mask2 = xAtt_mask;
    
    %x2(segBool) = nan;
    [~,fL] = min(abs(4.8e6- f));
    [~,fU] = min(abs(5.7e6 - f));

    for iF= fL:fU      
        
        %results are dB , dB/cm, [], y(iF) = a(iF)+b(iF)*2*x
    
        Flog = 10*log10(spectCorr(:,:,iF));
        FlogFilt = 10*log10(spectCorrFilt(:,:,iF));

        %Fcorr = 10*log10(Fcorr(iF,iXStart:iXEnd));
        %FcorrFilt = 10*log10(F_mask(iF,iXStart2:iXEnd2));

        options1 = [2*0.5*(endDepth + slDist) 0 0];
        
        xFit1 = 2*xAtt2./100;
        xFit2 = 2*xAtt_mask2./100;
        yFit1 = Flog;
        yFit2 = FlogFilt;
        
        coeffs = polyfit(xFit1(:),yFit1(:),1);
        a(iF) = coeffs(2);
        b(iF) = coeffs(1);


        validIdx  = ~isnan(xFit2) & ~isnan(yFit2);

        coeffs = polyfit(xFit2(validIdx),yFit2(validIdx),1);
        a2(iF) = coeffs(2);
        b2(iF) = coeffs(1);

        %[a2(iF),b2(iF),c2(iF),e2(iF,:),re2(iF)] = lsqn(xFit2,yFit2,options1);
        

        if 0%iF ==  21
            figure 
            subplot(1,2,1)
            plot(xFit1(:),yFit1(:),'k.','MarkerFaceColor','k')
            hold on 
            plot(xFit1(:),a(iF)+b(iF).*xFit1(:),'k-.','LineWidth',2)
            plot(xFit2(:),yFit2(:),'r.')
            plot(xFit2(validIdx),a2(iF)+b2(iF).*xFit2(validIdx),'r-*','LineWidth',2)
     
        end

    end

    b  = -b; 
    b2 = -b2;

    [~,iF0 ] = min(abs(f - 5.2e6))

    %results are dB/cm, dB/cm/MHz
    [a0,  alpha ,corr] = lsqn(f(fL:fU), b(fL:fU),[f(iF0) 5e5 2]); 
    [a02, alpha2,corr2] = lsqn(f(fL:fU), b2(fL:fU),[f(iF0) 5e5 2]); 
    
    if 0
        subplot(1,2,2)
        plot(f(fL:fU),b(fL:fU),'ko','MarkerFaceColor','k')
        hold on
        plot(f(fL:fU), a0+ alpha.*f(fL:fU),'k-.')
        plot(f(fL:fU), b2(fL:fU),'ro','MarkerFaceColor','r')
        plot(f(fL:fU), a02 + alpha2*f(fL:fU),'r-.')
    end

    resid = a0+ alpha.*f(fL:fU) - b(fL:fU);
    resid_Filt = a02 + alpha2.*f(fL:fU) - b2(fL:fU);

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