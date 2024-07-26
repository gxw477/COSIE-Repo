
clear 
close all 

%Image Analysis

topDirMaster = uigetdir;

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
kLength_COH = round(vsxParams.Receive(1).samplesPerWave);



  
%% 


fNames = ls([topDirMaster,'\BFimgData*']);
nImages =  size(fNames,1);

clearvars fNames

kIdxs = cell(nImages,1);

zSelect = input('Select depth of interest: ')/1e3;
[~, zIdx] = min(abs(rVals - zSelect));
axIdxs = zIdx-round(kLength_BSC_samples/2) : zIdx + round(kLength_BSC_samples/2) -1 ;


wOption = input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');

if wOption == 1 
    win = [0,ones(1,kLength_BSC_samples-2),0];
    wName = 'Rect';
elseif wOption == 2
    win = tukeywin(kLength_BSC_samples,0.1)';
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welsh';
end

saveDir = [topDirMaster,'/',wName,'/Z',num2str(round(zSelect*1e3)),'/'];

if ~exist(saveDir)
    mkdir(saveDir)
end


fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;

fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);
 
attSpeckle   = [0.524 , 0.09].*vsxParams.Trans.frequency;
%attWater = 0.00217*vsxParams.Trans.frequency^2; 

kCount = 1; 
rayCount= 0; 

imageIdxs = 1:20;

idxBool = 0.5*zSelect./5e-3 + (0:4:12);
idxBool = idxBool(idxBool<=20)

imageIdxs = imageIdxs(~ismember(1:nImages,idxBool))

nImages = length(imageIdxs);

for iImage = 1:nImages
    
    iImage

    load([topDirMaster,'\BFimgData',num2str(imageIdxs(iImage)),'.mat'])
     
    bM = bmode(iq',60);
    [~ , edge] = (max(bM,[],1));
    edgeYVal =  yVals(round(mean(edge)))

    attSpeckle_DB = 2*(zSelect*100 - edgeYVal*100)*attSpeckle(1);
    attComp_Speckle= 10^(attSpeckle_DB/10);
    

    imagesc((1:128).*lWidthF,yVals,bM); 
    colormap gray; 

    pause(0.5)

    lRl = 17; %input('Left ray line: ');
    rRl = 111;%input('Right ray line: ');

    kIdxs{iImage} = lRl:rRl;
   
    if length(axIdxs) ~= kLength_BSC_samples
        error('check bsc kernel length')
    end

    %% BSC + Coherence Calc

    bscLines = fullIM(kIdxs{iImage},axIdxs);
    
    if wOption == 1 || wOption == 2
        winMatrix = ones(size(bscLines)).*win;
        bscLines = bscLines.*winMatrix;
        spect = (fft(bscLines,[],2)).^2;
        
    elseif wOption == 3
        h = spectrum.welch;                  % Create a Welch spectral estimator.
        welchObj = psd(h,bscLines','Fs',fs);
        spect = (welchObj.Data').^2; % transpose because spectAveraging takes the vector the other way
    end
    
    spect2 = spect.*attComp_Speckle;
    spectAv = spectAveraging(spect2,kWidth_BSC_lines,oLap);

    nKsInAverage = size(spectAv,1);

    avSpectAll(kCount:(kCount+nKsInAverage-1),:) = spectAv;

    
    for iLine2 = 1:length(kIdxs{iImage})
        
       cohAll(iLine2 + rayCount ,:) = CoherenceAnalysisFN(squeeze(channelStack(kIdxs{iImage}(iLine2),axIdxs,:)));
       spectAll(iLine2 + rayCount, :) = spect2(iLine2,:);
       envAll(iLine2 + rayCount,:) = abs(envelope(bscLines(iLine2,:)));

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

input('check save folder ')

saveDir2 = [saveDir,'/COSIEoutput_adaptive'];

if ~exist(saveDir2)
    mkdir(saveDir2)
end


for sumIdx = 1:size(cohAll,2)

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



