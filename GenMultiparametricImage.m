


clear 
close all

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QA\';

%load verasonics param's
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

iImage =  input('Which image ? : ');

bfImgData = load([testDir,'BFimgData',num2str(iImage),'.mat']);


kLength_BSC_samples = 120;

fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
fVals = (0:kLength_BSC_samples-1).*fs/(kLength_BSC_samples);
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);

wOption = 2;%input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');

if wOption == 1 
    wName = 'Rect';
elseif wOption == 2
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welch';
end

adaptBool = 1;%input('Adaptive grid sizing? : ');

if adaptBool 
    adaptStr = '_adaptive';
else
    adaptStr = '';
end


depthVals = 15:5:40;
%regions of image with full aperture in effect
rayIdxs = 17:111;

sumIdx = 33;
kLength_COH_samples = 120;
kLength_COH_samples_Length = 0.5*kLength_COH_samples*vsxParams.sPerWaveInit*vsxParams.lambda;

EMLidx = 1;
kWidth = 5; 
oLap = 0.8;

COSIEall = cell(1,length(depthVals));
powAll = zeros(length(depthVals), length(rayIdxs));
cohAll = zeros(length(depthVals), length(rayIdxs));
segBoolAll = zeros(length(depthVals), length(rayIdxs));

for iDepth = 1:length(depthVals)

    depthSelect = depthVals(iDepth);

    dataDir = [testDir,wName,'\Z',num2str(depthSelect),'\'];
    % loop this over depths
    cInput = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);
    
    %centre frequency of test data
    powf0 = abs(cInput.spectAll(:,nF));
    powAll(iDepth,:) = powf0;

    speckleDir2 = [speckleDir,wName,'\Z',num2str(depthSelect),'\'];  
    speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(kLength_COH_samples),'\COSIEoutput',num2str(sumIdx),'.mat']);
    
    
    RMatTest = zeros(length(rayIdxs),sumIdx);
    
    [~ , cohTestKernelIdx] = min(abs(depthSelect*1e-3 - bfImgData.yVals));
    axIdxsCOH = (cohTestKernelIdx- kLength_COH_samples/2):(cohTestKernelIdx+kLength_COH_samples/2-1);

    for iLine = 1:length(rayIdxs)
       RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(rayIdxs(iLine),axIdxsCOH,:)));
    end

    cohTest = sum(RMatTest(:,1:sumIdx),2);
    cohAll(iDepth,:) = cohTest;

    segBool = cohTest > speckleCOSIE.redEML(1,EMLidx) & cohTest < speckleCOSIE.redEML(2,EMLidx);
  
    segBool1_cluster = ismember(rayIdxs, unique(cell2mat(idxClustering(rayIdxs2(segBool1),kWidth,oLap))))';

    segBoolAll(iDepth,:) = segBool1_cluster;

    COSIEall{iDepth} = speckleCOSIE;

end



samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample  + 1;
yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;

lambda = 1540/vsxParams.Trans.frequency*1e-6;
lWidth = lambda*(vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1)) ;
xVals = (1:128).*lWidth;


MultipleParametricCohImage(xVals,yVals,rayIdxs,bfImgData,depthVals,powAll,segBoolAll,kLength_BSC_samples,kLength_COH_samples)
