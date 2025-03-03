
clear 
close all 

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QAPht2\';
%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\G218L74_1\';

%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_HV_HTGC\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_NHV_NTGC\';

planeDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Pref\';

%load verasonics param's2
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);


sumIdx = 33;
iImage = input('Which Image ? : ');
%load test data
bfImgData = load([testDir,'BFimgData',num2str(iImage),'.mat']);

depthSelect = input('Depth of interest (mm) : ');
    
adaptBool = 1;%input('Adaptive grid sizing? : ');

if adaptBool 
    adaptStr = '_adaptive';
else
    adaptStr = '';
end

wOption = input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch \n 4 for Hanning: \n ');

if wOption == 1 
    wName = 'Rect';
elseif wOption == 2
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welch';
elseif wOption == 4
    wName = 'Hann';
end

saveDir_SEG = [testDir,'\',wName,'\'];

cohKlength = 120;
cohKlength_Length = 0.5*cohKlength*vsxParams.sPerWaveInit*vsxParams.lambda;

qaSNRdata = load([testDir,wName,'\Z',num2str(depthSelect),'\EnvStats',num2str(iImage),'.mat']);


samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample  + 1;
yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;

lambda = 1540/vsxParams.Trans.frequency*1e-6;
lWidth = lambda*(vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1)) ;
xVals = (1:128).*lWidth;

%Calculate kernel idxs for analysis
[~, depthIdx] = min(abs(yVals - depthSelect*1e-3));
kLength_BSC_samples = 120;
axIdxs = depthIdx-round(kLength_BSC_samples/2) : depthIdx + round(kLength_BSC_samples/2) -1 ;

%test data interface reflectivity
bm = bmode(bfImgData.iq',80);  
[~ , edge] = (max(bm,[],1));
edgeYVal = 0;%bfImgData.yVals(round(mean(edge)));

%spectral calculations
fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);

if ~exist([testDir,'/data_edgeCorr.mat'])
    
    %Calculate edge correction factor
    R0_test = planarReflectorEstimates(testDir,planeDir,0);
    Ttest = (1-mean(R0_test));
    
    R0_speckle= planarReflectorEstimates(speckleDir,planeDir,1);
    Tspeckle = (1-mean(R0_speckle));

    edgecorr = (Tspeckle/Ttest)^2;
    save([testDir,'/data_edgeCorr.mat'],'edgecorr')
    
else
    load([testDir,'/data_edgeCorr.mat'])
    %
end
 
input('Check which attenuation value you are using : ')

%QA phantom
attTest     = [0.579 , 0.955].*vsxParams.Trans.frequency;

%Emma Liver 
%attTest     = [0.562 , 0.8].*vsxParams.Trans.frequency;


attSpeckle  = [0.524 , 0.09].*vsxParams.Trans.frequency;
%attP        = [0.53  , 0.003].*vsxParams.Trans.frequency;


%ground truth BSC values: reference 
mu0_ref = 3.86e-4/(3^3.5); %(cm^{-1}sr^{-1}) 
bscSpeckleBf_ref = mu0_ref * 100 * vsxParams.Trans.frequency^3.5 ; %(m^{-1}sr^{-1})
bscSpeckleSTD_ref = 0.2*bscSpeckleBf_ref;

%ground truth BSC values: test
mu0_test = 3.25e-4/(3^3.6);%(cm^{-1}sr^{-1}) 
bscSpeckleBf_test = mu0_test * 100 * vsxParams.Trans.frequency^3.6; %(m^{-1}sr^{-1})
bscSpeckleSTD_test = 0.2*bscSpeckleBf_test;



%regions of image with full aperture in effect
xBool = 17:112;

%calculate coherence properties
[~ , cohTestKernelIdx] = min(abs(depthSelect*1e-3 - bfImgData.yVals));

RMatTest = zeros(length(xBool),sumIdx);
axIdxsCOH = (cohTestKernelIdx- cohKlength/2):(cohTestKernelIdx+cohKlength/2-1);

for iLine = 1:length(xBool)
   RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(xBool(iLine),axIdxsCOH,:)));
end

cohTest = sum(RMatTest(:,1:sumIdx),2);

snrTest = qaSNRdata.envMean./qaSNRdata.envStd;

%% Generate parametric images using first 3 segmentation points on EML


rayIdxs = (1:128);    
rayIdxs2 = rayIdxs(xBool);  
EMLidx = 3;
kWidth = 5; 
oLap = 0.8;

allDepths = 15:5:50;

segBoolBIG = zeros(length(rayIdxs2),length(allDepths));
powf0_BIG = segBoolBIG;

axIdxs_BIG = [];

for iDepths = 1:length(allDepths)

    dzT = allDepths(iDepths)*0.1 - edgeYVal*100; %cm 
    attTest_DB = 2*(dzT)*attTest(1);
    attComp_Test =   10^(attTest_DB/10);

    [~ , cohTestKernelIdx] = min(abs(allDepths(iDepths)*1e-3 - bfImgData.yVals));

    RMatTest = zeros(length(xBool),sumIdx);
    axIdxsCOH = (cohTestKernelIdx- cohKlength/2):(cohTestKernelIdx+cohKlength/2-1);
    
    for iLine = 1:length(xBool)
       RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(xBool(iLine),axIdxsCOH,:)));
    end
    
    cohTest = sum(RMatTest(:,1:sumIdx),2);

    speckleDir2 = [speckleDir,wName,'\Z',num2str(allDepths(iDepths)),'\'];  
    dataDir = [testDir,wName,'\Z',num2str(allDepths(iDepths)),'\'];
    cInput = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);

    %load COSIE data
    speckleSNRdata= load([speckleDir2 , 'envData_COSIE' ]) ;

    speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(cohKlength),'\COSIEoutput',num2str(sumIdx),'.mat']);
    cInput_Depth = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);

    powf0 = abs(cInput_Depth.spectAll(:,nF));

    powerSeg      = COVsegmentation_sK(cohTest,speckleCOSIE.EML,(powf0),kWidth,oLap);

    %mean speckle power
    specklePOWER_MEAN = mean(speckleCOSIE.powf0);
    
    bscEstimate = (powf0./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;
    
    segBool1 = cohTest > speckleCOSIE.redEML(1,EMLidx) & cohTest < speckleCOSIE.redEML(2,EMLidx);    
    segBool1_cluster = ismember(rayIdxs2, unique(cell2mat(idxClustering(rayIdxs2(segBool1),kWidth,oLap))))';
    
    powf0_BIG(:,iDepths) = bscEstimate;

    segBoolBIG(:,iDepths) = segBool1_cluster;
    
    axIdxs_BIG = [axIdxs_BIG,axIdxsCOH];
        
end

%axIdxs_BIG = fIdx(1):fIdx(2);


[ax1,ax2,cB ] = bmCohCOSIEparImage_mDepth(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs_BIG,kLength_BSC_samples,powf0_BIG,segBoolBIG,xBool,speckleCOSIE.pctSeg2(EMLidx),'Coherence');
xlim([-15 15])
ax1.XTick = [-15 :5:15];
ax1.YTick = [5:10:55];
ax1.YLim = [5 55];
set(ax1,'FontSize',14)
set(ax1,'FontWeight','normal')
saveas(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage)])
saveas(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage),'.jpg'])
savefigPDF_Crop(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage)])


%bmCohCOSIEparImage_mDepth(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs_BIG,kLength_BSC_samples,powf0_BIG,segBoolBIG,xBool,speckleCOSIE.pctSeg2(EMLidx),'SNR')


%save([saveDir_SEG,'\SegResults.mat'],'bscSpeckleBf_ref','bscSpeckleSTD_ref','bscEstimate_COH_COSIE',...
%    'bscEstimate_ENV_COSE','bscEstimate_ENV_WEIGHT','bscEstimate_COH_WEIGHT')


