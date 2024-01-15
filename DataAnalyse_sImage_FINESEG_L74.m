
for  sumIdx = 1:33
    for iImage = 1:3
    
        clearvars -except iImage sumData sumIdx
        close all 
        
        %topDir = [uigetdir(cd,'Select Analysis directory'),'\'];
        topDir = ['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74_QuadInterp\QASpeckle2\\'];
        
        adaptBool = 1;%input('Adaptive grid sizing? : ');
        
        if adaptBool 
            adaptStr = '_adaptive';
        else
            adaptStr = '';
        end
        
        speckleDir = ['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74_QuadInterp\Speckle\'];
        
        vsxParams = load([topDir,'\VSXoutput.mat']);
        vsxParams2 = load([topDir,'\VSXoutput.mat']);
        
        isequal(vsxParams.TGC,vsxParams2.TGC)
        
        %iImage = 3;%input('Image # : ');
        
        load([topDir,'COSIEinput',num2str(iImage),'.mat'])
        load([topDir,'BFimgData',num2str(iImage),'.mat'])
        
        
        %% 1. Set up variables 
        
        samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample  + 1;
        
        yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;
        
        [~, bscFocIdx] = min(abs(yVals - vsxParams.TX(1).focus*vsxParams.lambda));
        kLength_BSC_samples = 120;
        axIdxs = bscFocIdx-round(kLength_BSC_samples/2) : bscFocIdx + round(kLength_BSC_samples/2) -1 ;
        
        fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
        fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
        df = fVals(2)-fVals(1);
        nF = round(vsxParams.Trans.frequency*1e6/df);
         
        powf0 = abs(spectAll(:,nF));
        
        %sumIdx = size(cohAll,2);
        
        
        %% 2. Compute segmentation based on EML
       
        imgDir = [topDir,'/SegResults_',num2str(iImage),adaptStr,'/SumIdx_',num2str(sumIdx)];
        
        if ~exist(imgDir)
            mkdir(imgDir)
        end
        
        speckleCOSIE = load([speckleDir,'\COSIEoutput',adaptStr,'\COSIEoutput',num2str(sumIdx),'.mat']);
        
        
        cohSum = sum(cohAll(:,1:sumIdx),2);
        
        
        %This step is for interest only, just want to see how the EML's differ
        cosieParams.dTH= 0.1;
        cosieParams.APsize = sumIdx;
        [bscSurface,segSurface,EML,pctSeg1,redEML,pctSeg2] = COSIE(cohSum,powf0,cosieParams);
        
        
        transpData = zeros(size(fullIM));
        colorData = nan.*zeros(size(fullIM));
        
        kWidth = 5; 
        oLap = 0.8;
        
        tGCV = linspace(0.3,1,size(fullIM,2));
        tGCM = repmat(tGCV,size(fullIM,1),1);
        
        iSegPct = 1;
            
        %1 = keep  , 0 = segment
        segBool = cohSum > speckleCOSIE.EML(1,iSegPct) & cohSum < speckleCOSIE.EML(2,iSegPct);
        
        outSeg = (1-length(find(segBool))/length(segBool))*100;
        
        %currently a bit off, should be jumping where we can't make a full
        %kernel widthwise
        %spectOut = spectAveraging(spectAll(segBool), kWidth, oLap);
        
        compImage2 = log(abs(compImage'));
        
        compImage2(isinf(compImage2(:))) = 0;
        
        compImage2 = compImage2.*tGCM';
        
        
        %% 3. Plot BSC/b-mode image
        lWidth = (Trans.ElementPos(2,1)-Trans.ElementPos(1,1))*lambda;
        ax1 = axes;
        
        B80 = bmode(iq',80);
        
        imagesc(ax1,lWidth.*(1:128),yVals,B80);
        hold on 
        plot(ax1,lWidth.*[1 size(fullIM,1)],yVals(bscFocIdx).*[1 1],'-','color','red')
        plot(ax1,lWidth.*120.*[1 1], [yVals(axIdxs(1)) yVals(axIdxs(end))],'-','color','red')
        xlabel('Lateral Position (m)')
        axis equal
        ylabel('Axial Position (m)')
        
        ax2 = axes;
        powfLong = repmat(powf0,1,size(axIdxs,2));
        
        hideVectorLong = repmat(segBool,1,size(axIdxs,2));
        
        colorData(rayIdxs,axIdxs) = log(powfLong); 
        
        transpData(rayIdxs,axIdxs) = 0.6.*hideVectorLong; 
        
        imagesc(ax2 ,lWidth.*(1:128) , yVals, colorData','AlphaData',transpData')    
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
        
        title(ax1,['Seg = ',num2str(outSeg),'%'])
        axis tight 
        axis equal
        
        fname2 = [imgDir,'/ParametricImage'];
        savefig(fname2)
        saveas(gcf,fname2,'png')
        
        
        
        cohValues2 = sum(RMat(:,:,1:sumIdx),3);
        F = griddedInterpolant(xQcoh,yQcoh,cohValues2);
        CoherenceQ = F(xQbm,yQbm);
            
        
        figure 
        imagesc(lWidth.*(1:128) , yVals, CoherenceQ');
        xlabel('Lateral Position (m)')
        ylabel('Axial Position (m)')
        hold on 
        plot(lWidth.*[1 size(fullIM,1)],yVals(bscFocIdx).*[1 1],'-','color','red')
        plot(lWidth.*120.*[1 1], [yVals(axIdxs(1)) yVals(axIdxs(end))],'-','color','red')
        colormap gray
        axis tight 
        axis equal
        
        fname2 = [imgDir,'/CoherenceImage'];
        savefig(fname2)
        saveas(gcf,fname2,'png')
        
        %% 4. Plot beam features at the focus
        
        figure
        subplot(3,1,1)
        plot(rayIdxs,cohSum)
        hold on
        plot([1 128],[speckleCOSIE.EML(2,1) speckleCOSIE.EML(2,1)],'-','Color','red')
        plot([1 128],[speckleCOSIE.EML(1,1) speckleCOSIE.EML(1,1)],'-','Color','red')
        ylim([-1.5*max(speckleCOSIE.EML(2,1)) 1.5*max(speckleCOSIE.EML(2,1))])
        xlabel('Rayline')
        ylabel('Summed Coherence')
        xlim([1 128])
        
        
        subplot(3,1,2)
        
        
        focusEnv = abs(envelope(fullIM(:,axIdxs)'));
        mEnv = mean(focusEnv,1);
        sEnv = std(mEnv,1);
        plot(mEnv,'-.','Color','k')
        set(gca,'YScale','log')
        yticks(linspace(min(mEnv)*0.5,max(mEnv)*1.1,5));
        ylim([min(mEnv)*0.5 max(mEnv)*1.1])
        
        xlabel('Rayline')
        ylabel('Mean envelope Amplitude at focus')
        xlim([1 128])
        
        
        subplot(3,1,3)
        plot(rayIdxs(segBool),powf0(segBool),'x','Color','red')
        hold on 
        plot(rayIdxs(~segBool),powf0(~segBool),'o','Color','red')
        xlabel('Ray line')
        ylabel('Backscattered power (a.u.)')
        
        mPow = max(powf0);
        xlim([1 128])
        ylim([-0.1*mPow 1.1*mPow])
        
        set(gca,'YScale','log')
        
        
        fname2 = [imgDir,'\FocLines'];
        savefig(fname2)
        saveas(gcf,fname2,'png')
        
        
        %% 7. Repeat 5. for more seg pctgs 
        outSegPrev = 0; 
        
        
        kWidth = 5;
        oLap = 0.8;
        kIdxs = idxClustering(1:size(spectAll,1),kWidth,oLap);
        nKernels = size(kIdxs,2);
        bscValuesIK = zeros(1,nKernels);
        stdValuesIK = zeros(1,nKernels);
        
        nPossibleKernels = size(idxClustering(1:length(cohSum),kWidth,oLap),2);
        
        
        for iKernels = 1:nKernels
            bscValuesIK(iKernels) = abs(mean(spectAll(kIdxs{iKernels},nF)));
        end
        
        segPct =  100*(1- nKernels/nPossibleKernels);
            
        bscEstimate2(:,1) = [mean(bscValuesIK) ;std(bscValuesIK); 0];
        
        segPct = 0;
        segPctIdx = 2;
        
        while segPct < 50
        
            if segPct > 8 && segPct < 9
                segPct;
            end
        
            segBool = cohSum > speckleCOSIE.EML(1,segPctIdx) & cohSum < speckleCOSIE.EML(2,segPctIdx);
            
            segIdxs = find(segBool);
            
            kIdxs = idxClustering(segIdxs,kWidth,oLap);
        
            nKernels = size(kIdxs,2);
        
            bscValuesIK = zeros(1,nKernels);
            stdValuesIK = zeros(1,nKernels);
        
            for iKernels = 1:nKernels
                bscValuesIK(iKernels) = abs(mean(spectAll(kIdxs{iKernels},nF)));
            end
        
            segPct =  100*(1- nKernels/nPossibleKernels);
        
            bscEstimate2(:,segPctIdx) = [mean(bscValuesIK) ;std(bscValuesIK);segPct];
            
            segPctIdx = segPctIdx+1;
        
        end
        
        
        nValues = bscEstimate2(3,:) * length(cohSum);
        
        figure 
        subplot(1,3,1)
        errorbar(bscEstimate2(3,:),bscEstimate2(1,:),bscEstimate2(2,:),'-.')
        set(gca,'YScale','log')
        xlabel('Segmentation %')
        ylabel('Backscattered power (a.u.)')
        title('Power')
        xlim([-5 max(bscEstimate2(3,:))+10])
        yMax = max(bscEstimate2(1,:));
        yMin = min(bscEstimate2(1,:));
        ylim([yMin*0.5 yMax*4])
        tickVector = round(10.^linspace(log10(10),log10(4*yMax),10));
        yticks(tickVector)
        
        subplot(1,3,2)
        plot(bscEstimate2(3,:),bscEstimate2(2,:),'-x')
        set(gca,'YScale','log')
        xlabel('Segmentation %')
        ylabel('Standard Deviation of BSC (a.u.)')
        title('Standard Deviation')
        xlim([-5 max(bscEstimate2(3,:))+10])
        yMax = max(bscEstimate2(2,:));
        yMin = min(bscEstimate2(2,:));
        ylim([yMin*0.5 yMax*4])
        tickVector = round(10.^linspace(log10(10),log10(4*yMax),10));
        yticks(tickVector)
        
        
        subplot(1,3,3)
        COV = bscEstimate2(2,:)./bscEstimate2(1,:);
        plot(bscEstimate2(3,:),COV,'-x')
        set(gca,'YScale','log')
        xlabel('Segmentation %')
        ylabel('C.O.V.')
        title('Coefficient of Variation')
        xlim([-5 max(bscEstimate2(3,:))+10])
        xlim([-5 max(bscEstimate2(3,:))+10])
        yMax = max(COV);
        yMin = min(COV);
        ylim([yMin*0.5 yMax*2])
        
        tickVector = linspace(0,3,12);
        yticks(tickVector)
        
        fname2 = [imgDir,'\VariabilityEstimates'];
        savefig(fname2)
        saveas(gcf,fname2,'png')
        
    end
end