
clear 
close all 

%Image Analysis

transSwitch = input('0 for L74 \n1 for C1-6D \n : ');

if transSwitch == 0
    topDirMaster =  'C:\Users\gwest\Documents\MATLAB\ElastPhtL74\Img1-4Dir\';
elseif transSwitch == 1 
    topDirMaster =  'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE_StudyData\ElastPht\BFimgDataTGCCorr\';
end

vsxParams = load([topDirMaster,'/VSXoutput.mat']);


%% Define BSC kernel properties
lambda = vsxParams.lambda;

%dtheta = vsxParams.Angle(2)-vsxParams.Angle(1);

kWidth_BSC_lines = 5;
kLength_BSC_samples = 120;
oLap = .80;

lWidthF = lambda*(vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1)) ;

samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample + 1;

%already defined, but to remind you
rVals = lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;

[~, bscFocIdx] = min(abs(rVals - vsxParams.TX(1).focus*lambda));


%% Define coherence kernel properties 

%kernel length (wavels)
nWavels = 3;
%kernel length (samples)
kLength_COH = kLength_BSC_samples*(4/4);

%% 

fNames = ls([topDirMaster,'\BFimgData*']);
nImages =  size(fNames,1);

clearvars fNames

kIdxs = cell(nImages,1);



for zSelect = (15:5:50).*1e-3

    [~, zIdx] = min(abs(rVals - zSelect));
    axIdxsBSC = zIdx-round(kLength_BSC_samples/2) : zIdx + round(kLength_BSC_samples/2) -1 ;
    axIdxsCOH = zIdx-round(kLength_COH/2) : zIdx + round(kLength_COH/2) -1 ;
    
    
    
    wOption = 4;%input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');
    
    if wOption == 1 
        win = [0,ones(1,kLength_BSC_samples-2),0];
        wName = 'Rect';
    elseif wOption == 2
        win = tukeywin(kLength_BSC_samples,0.1)';
        wName = 'Tukey';
    elseif wOption == 3
        wName = 'Welsh';
    elseif wOption == 4 
        win = hann(kLength_BSC_samples)';
        wName = 'Hann';
    end
    
    saveDir = [topDirMaster,'/',wName,'/Z',num2str(round(zSelect*1e3)),'/'];
    
    if ~exist(saveDir)
        mkdir(saveDir)
    end
    
    
    fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
    
    fVals = (0:length(axIdxsBSC)-1).*fs/(length(axIdxsBSC));
    df = fVals(2)-fVals(1);
    nF = round(vsxParams.Trans.frequency*1e6/df);
     
    attSpeckle   = [0.524 , 0.09].*vsxParams.Trans.frequency;
    %attWater = 0.00217*vsxParams.Trans.frequency^2; 
    
    kCount = 1; 
    rayCount= 0; 
    
    if transSwitch == 0
        %all possible image idxs
        imageIdxsAll = 1:nImages;
        
        %number of beam translations IN FOLDER
        nTransl = 4%input('Number of beam translations in folder : ');
        %number of frames/Sets 
        nFrames = nImages/nTransl;%input('Number of Repeats : ');
        %increment between images
        incZ = 5e-3;
        %image idxs within frame
        imageIdxsFrame = 1:nTransl;
        
        if round(nImages/nTransl)~= nImages/nTransl
            error('wrong number of images')
        end
        
        
        %reject if the analysis kernel is on the water or on the edge 
        edgeBoolFile = imageIdxsFrame < zSelect./incZ ;
        %reject if the analysis kernel includes a reverb
        revbBoolFile = 0.5.*zSelect./incZ ~= imageIdxsFrame;
        
        %which numbers in 1:nTransl are neither half or equal to the selected depth
        %or in the water path
        boolFile = edgeBoolFile & revbBoolFile ;
        imageIdxsFrameKeep = imageIdxsFrame(boolFile);
        %duplicate to get enough rows for all frames
        imageIdxsFrameKeep2 = repmat(imageIdxsFrameKeep,nFrames,1);

        %add iFrame*nTransl to each row to get all available images
        vec1 = nTransl.*[0:(nFrames-1)]';
        mat1 = repmat(vec1,[1,size(imageIdxsFrameKeep2,2)]);
        
        imageIdxsFrameKeep3 = imageIdxsFrameKeep2 + mat1
        %Apply across nFrames to find all the indices
        imageIdxsAll = sort(imageIdxsFrameKeep3(:));
    
    elseif transSwitch == 1
        
        %we're using all of them because they're recorded at the same
        %distance from the phantom
        imageIdxsAll = 1:nImages;

    end
    
    nImages2 = length(imageIdxsAll);
    
    for iImage = 1:nImages2
        
        iImage
    
        load([topDirMaster,'\BFimgData',num2str(imageIdxsAll(iImage)),'.mat'])
         
        bM = bmode(iq',60);
        
        [~,edgeIdx1] = min(abs(yVals-3e-3));
        [~,edgeIdx2] = min(abs(yVals-3.6e-2));

        [~ , edgeIdxActual] = max(bM(1:edgeIdx2,64));
        edgeYVal =  yVals(edgeIdxActual+edgeIdx1-1);

        figure 
        imagesc(1:128,yVals,bM)
        colormap gray 
        title(['Edge at ',num2str(edgeYVal*1e3)])

        close all
    
        attSpeckle_DB = 2*(zSelect*100 - edgeYVal*100)*attSpeckle(1);
        attComp_Speckle= 10^(attSpeckle_DB/10);

        temp(iImage) = edgeYVal;
        
        if transSwitch == 0
            imagesc((1:128).*lWidthF,yVals,bM); 
            colormap gray; 
    
            pause(0.5)
    
            lRl = 17; %input('Left ray line: ');
            rRl = 111;%input('Right ray line: ');
    
        elseif transSwitch == 1
            
            angleBool = abs(vsxParams.Angle) <= 0.1;
            angleBoolIdx = find( angleBool);
            lRl = angleBoolIdx(1);
            rRl = angleBoolIdx(end);
            
        end
        
        kIdxs{iImage} = lRl:rRl;
       
        if length(axIdxsBSC) ~= kLength_BSC_samples
            error('check bsc kernel length')
        end
    
        %% BSC + Coherence Calc
    
        bscLines1 = fullIM(kIdxs{iImage},axIdxsBSC);
        
        if wOption == 1 || wOption == 2 || wOption == 4
            winMatrix = ones(size(bscLines1)).*win;
            bscLines2 = bscLines1.*winMatrix;
            spect = (fft(bscLines2,[],2)).^2;
            
        elseif wOption == 3
            h = spectrum.welch;                  % Create a Welch spectral estimator.
            welchObj = psd(h,bscLines1','Fs',fs);
            spect = (welchObj.Data').^2; % transpose because spectAveraging takes the vector the other way
        end
        
        spect2 = spect.*attComp_Speckle;
        spectAv = spectAveraging(spect2,kWidth_BSC_lines,oLap);
    
        nKsInAverage = size(spectAv,1);
    
        avSpectAll(kCount:(kCount+nKsInAverage-1),:) = spectAv;
    
        
        for iLine2 = 1:length(kIdxs{iImage})
            
           cohAll(iLine2 + rayCount ,:) = CoherenceAnalysisFN(squeeze(channelStack(kIdxs{iImage}(iLine2),axIdxsCOH,:)));
           spectAll(iLine2 + rayCount, :) = spect2(iLine2,:);
           envAll(iLine2 + rayCount,:) = abs(envelope(bscLines1(iLine2,:)));
    
        end
        
        kCount = kCount + nKsInAverage; 
        rayCount = rayCount + length(kIdxs{iImage});
    
    end
    
    
    powf0 = abs(spectAll(:,nF));
    
    mEnv = mean(envAll,2);
    sEnv = std(envAll,0,2);
    snr = mEnv./sEnv;
    
    [bscSurface,segSurface,EML,pctSeg1,redEML,pctSeg2] = COSIE_adaptiveGrid(snr,powf0);
    
    save([saveDir,'envData_COSIE.mat'],'snr','EML','bscSurface','pctSeg1','pctSeg2','redEML','segSurface')
    
    avCohAll = spectAveraging(cohAll,kWidth_BSC_lines,oLap);
    
    
    cohMubCorr = zeros(1,size(cohAll,2));
    
    saveDir2 = [saveDir,'/COSIEoutput_adaptive',num2str(kLength_COH)];
    
    if ~exist(saveDir2)
        mkdir(saveDir2)
    end
    
   
    for sumIdx = [16 ,size(cohAll,2)]
    
        sumIdx/size(cohAll,2)
    
        cohSum = sum(cohAll(:,1:sumIdx),2);
    
        cosieParams.APsize = sumIdx;
        cosieParams.dTH = sumIdx/0.5e3;
        cosieParams.GT = (mean(powf0));
        
        [bscSurface,segSurface,EML,pctSeg1,redEML,pctSeg2] = COSIE_adaptiveGrid(cohSum,powf0);
        
        thVector = sort(cohSum);
    
        save([saveDir2,'/COSIEoutput',num2str(sumIdx),'.mat'],'EML','bscSurface','pctSeg1','pctSeg2','redEML','segSurface','thVector','powf0','cohSum')
    
        %plot(EML(1,:),EML(2,:),'.')
        %close all
        
    
        R = corrcoef(cohSum,abs(powf0));
    
        cohMubCorr(sumIdx) = R(2);
    
    end

end

