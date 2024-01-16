
clear 
close all 

%Image Analysis

topDir = 'PhantomExperimentsL74_QuadInterp/Speckle/';

vsxParams = load([topDir,'/VSXoutput.mat']);


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
kLength_COH = round(vsxParams.Receive(1).samplesPerWave);



  
%% 
nImages =  12;
kIdxs = cell(nImages,1);

axIdxs = bscFocIdx-round(kLength_BSC_samples/2) : bscFocIdx + round(kLength_BSC_samples/2) -1 ;

fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;

fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);
 

kCount = 1; 
rayCount= 0; 

rayIdxFname = [topDir,'/rayIdxs.mat'];

for iImage = 1:nImages

    iImage

    load([topDir,'BFimgData',num2str(iImage),'.mat'])
    
    if 1 %~exist([rayIdxFname])
        
        bM = bmode(iq,30);
    
        imagesc(bM'); 
        colormap gray; 

        pause(0.5)

        lRl = 17;%input('Left ray line: ');
        rRl = 111;%input('Right ray line: ');

        kIdxs{iImage} = lRl:rRl;
    
    else
        
        load(rayIdxFname)
        
    end


   
    if length(axIdxs) ~= kLength_BSC_samples
        error('check bsc kernel length')
    end

    %% BSC + Coherence Calc

    bscLines = fullIM(kIdxs{iImage},axIdxs);

    spect = fft(bscLines,[],2);

    
    spectAv = spectAveraging(spect,kWidth_BSC_lines,oLap);

    nKsInAverage = size(spectAv,1);

    avSpectAll(kCount:(kCount+nKsInAverage-1),:) = spectAv;

    
    for iLine2 = 1:length(kIdxs{iImage})
        
       cohAll(iLine2 + rayCount ,:) = CoherenceAnalysisFN_new(squeeze(channelStack(kIdxs{iImage}(iLine2),axIdxs,:)));
       spectAll(iLine2 + rayCount, :) = spect(iLine2,:);
       envAll(iLine2 + rayCount,:) = abs(envelope(bscLines(iLine2,:)));

    end
    
    kCount = kCount + nKsInAverage; 
    rayCount = rayCount + length(kIdxs{iImage});

end


save([topDir,'envData.mat'],'envAll')

avCohAll = spectAveraging(cohAll,kWidth_BSC_lines,oLap);

powf0 = abs(spectAll(:,nF));

cohMubCorr = zeros(1,size(cohAll,2));

for sumIdx = 1:size(cohAll,2)

    cohSum = sum(cohAll(:,1:sumIdx),2);

    cosieParams.APsize = sumIdx;
    cosieParams.dTH = sumIdx/0.5e3;
    cosieParams.GT = (mean(powf0));
    
    [bscSurface,segSurface,EML,pctSeg1,redEML,pctSeg2,params] = COSIE_adaptiveGrid(cohSum,powf0,cosieParams);
    
    thVector = sort(cohSum);

    save([topDir,'COSIEoutput_adaptive/COSIEoutput',num2str(sumIdx),'.mat'],'EML','bscSurface','pctSeg1','pctSeg2','redEML','segSurface','params','thVector')

    R = corrcoef(cohSum,abs(powf0));

    cohMubCorr(sumIdx) = R(2);

end
