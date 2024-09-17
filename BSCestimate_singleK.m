
clear 
close all 

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QAPht1\';
%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\G218L74_1\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver_NHV_NTGC\';


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

wOption = 2;%input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');

if wOption == 1 
    wName = 'Rect';
elseif wOption == 2
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welch';
end

cohKlength = 120;
cohKlength_Length = 0.5*cohKlength*vsxParams.sPerWaveInit*vsxParams.lambda;

speckleDir2 = [speckleDir,wName,'\Z',num2str(depthSelect),'\'];  
%load COSIE data
speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(cohKlength),'\COSIEoutput',num2str(sumIdx),'.mat']);
speckleSNRdata= load([speckleDir2 , 'envData_COSIE' ]) ;


dataDir = [testDir,wName,'\Z',num2str(depthSelect),'\'];
cInput = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);
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

%centre frequency of test data
powf0 = abs(cInput.spectAll(:,nF));


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
 
input('Check which attenuation value you are using')

%QA phantom
%attTest     = [0.579 , 0.955].*vsxParams.Trans.frequency;

%Emma Liver 
attTest     = [0.562 , 0.8].*vsxParams.Trans.frequency;


attSpeckle  = [0.524 , 0.09].*vsxParams.Trans.frequency;
%attP        = [0.53  , 0.003].*vsxParams.Trans.frequency;



dzT = depthSelect*0.1 - edgeYVal*100;
attTest_DB = 2*(dzT)*attTest(1);
attComp_Test =   10^(attTest_DB/10);

%mean speckle power
specklePOWER_MEAN = mean(speckleCOSIE.powf0);
speckle_STD  = std(speckleCOSIE.powf0);

%mean test power
testPOWER_MEAN = mean(powf0*attComp_Test);
testPOWER_STD = std(powf0*attComp_Test);

%ground truth BSC values: reference 
mu0_ref = 3.86e-4/(3^3.5); %(cm^{-1}sr^{-1}) 
bscSpeckleBf_ref = mu0_ref * 100 * vsxParams.Trans.frequency^3.5 ; %(m^{-1}sr^{-1})
bscSpeckleSTD_ref = 0.2*bscSpeckleBf_ref;

%ground truth BSC values: test
mu0_test = 3.25e-4/(3^3.6);%(cm^{-1}sr^{-1}) 
bscSpeckleBf_test = mu0_test * 100 * vsxParams.Trans.frequency^3.6; %(m^{-1}sr^{-1})
bscSpeckleSTD_test = 0.2*bscSpeckleBf_test;



%regions of image with full aperture in effect
xBool = 17:111;

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

saveDir_SEG = [testDir,wName,'\Z',num2str(depthSelect),'\SegResults_',num2str(iImage),adaptStr,num2str(cohKlength),'\SumIdx_',num2str(sumIdx),'\'];
  
if ~exist(saveDir_SEG,'dir')
    mkdir(saveDir_SEG)
end



EMLidx = 1;
kWidth = 5; 
oLap = 0.8;

segBool1 = cohTest > speckleCOSIE.redEML(1,EMLidx) & cohTest < speckleCOSIE.redEML(2,EMLidx);

%bmCohCOSIEparImage takes segBool not indices as an argument, so we'll
%convert to idx and back to bool to give a segbool that reflects the actual
%useable lines. We want the second output argument of idxsClustering to 
%just get the unique lines instead of a cell array (output 1).

%All lines
rayIdxs = (1:128);

%Lines in centre of image
rayIdxs2 = rayIdxs(xBool);


segBool1_cluster = ismember(rayIdxs2, unique(cell2mat(idxClustering(rayIdxs2(segBool1),kWidth,oLap))))';


bmCohCOSIEparImage(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs,axIdxsCOH,powf0,segBool1_cluster,xBool,speckleCOSIE.pctSeg2(EMLidx),'Coherence')
saveas(gcf,[saveDir_SEG,'ParametricImage_COH'])


segBool2 = snrTest > speckleSNRdata.redEML(1,EMLidx) & snrTest < speckleSNRdata.redEML(2,EMLidx);
bmCohCOSIEparImage(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs,axIdxsCOH,powf0,segBool2,xBool,speckleSNRdata.pctSeg2(EMLidx),'SNR')
saveas(gcf,[saveDir_SEG,'ParametricImage_SNR'])

depthFeaturePlot(speckleCOSIE,speckleSNRdata,EMLidx,lWidth,cohTest,snrTest,20)
saveas(gcf,[saveDir_SEG,'DepthImage'])





%%

nEMLpoints = size(speckleCOSIE.EML,2);
%bscCOSIE = zeros(2,nEMLpoints+1);

powerSeg      = COVsegmentation_sK(cohTest,speckleCOSIE.EML,(powf0),kWidth,oLap);


bscEstimate = zeros(size(powerSeg,2),4);
bscEstimate(:,3) = powerSeg(3,:);
bscEstimate(:,4) = powerSeg(4,:);

bscEstimate(:,1) = (powerSeg(1,:)./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;
bscEstimate(:,2) = ((powerSeg(2,:))./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;

errs = zeros(size(bscEstimate,1),6);
%[perrSi,perrSp,perralphaT,perralphaR,perrbscR , TOTAL] x nSegmentations

for i = 1:size(bscEstimate,1)

    errs(i,:) = errorPropagationBSC(powerSeg(1,i),powerSeg(2,i),specklePOWER_MEAN,speckle_STD,dzT,attTest(2),2e-2,attSpeckle(2),bscSpeckleBf_ref,0.2*bscSpeckleBf_ref,bscEstimate(i,1));
    %bscEstimate(i,3) = errs(i,6);

end



orderSlices = [3,4,2,5,1];

legendCell = cell(5,1);
legendCell(1) = {'\epsilon({S_i})'};
legendCell(2) = {'\epsilon(S_p)'}; 
legendCell(3) = {'\epsilon(\alpha_{T})'};
legendCell(4) = {'\epsilon(\alpha_{R})'};
legendCell(5) = {'\epsilon(\mu_{R})'};


figure 
pie(errs(1,orderSlices)./sum(errs(1,1:5)),[0 5 0 5 0])
%legend({'\epsilon(\alpha_{T})','\epsilon(\alpha_R)','\epsilon(S_p)','\epsilon(\mu_{T})','\epsilon({S_i})'})
legend(legendCell(orderSlices))


%savefig(gcf,['BSC_EstimateFigure/S_',num2str(sumIdx)])
%saveas(gcf,['BSC_EstimateFigure/S_',num2str(sumIdx),'.jpg'])


%% more parametrics? 

%contBool = input('More maps? ');

if 0
    
    EMLvideoMaker([saveDir_SEG,'\video\'],55,cohTest,speckleCOSIE,lWidth,yVals,bfImgData,depthIdx,axIdxs,powf0,xBool,sumIdx)

end


%% Calculations for segmentations vs BSC estimate with both methods and both input vbls

%close all

testEnvelope = load([dataDir,'\EnvStats',num2str(iImage),'.mat']);
    


powerSeg_COH_COSIE = COVsegmentation_sK(cohTest,speckleCOSIE.redEML,powf0,kWidth,oLap);
bscEstimate_COH_COSIE = zeros(size(powerSeg_COH_COSIE,2),4);
bscEstimate_COH_COSIE(:,3) = powerSeg_COH_COSIE(3,:);
bscEstimate_COH_COSIE(:,4) = powerSeg_COH_COSIE(4,:);
bscEstimate_COH_COSIE(:,1) = (powerSeg_COH_COSIE(1,:)./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;
bscEstimate_COH_COSIE(:,2) = ((powerSeg_COH_COSIE(2,:))./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;


snrEnv = testEnvelope.envMean./testEnvelope.envStd;
powerSeg_ENV_COSIE = COVsegmentation_sK(snrEnv,speckleSNRdata.redEML,powf0,kWidth,oLap);
bscEstimate_ENV_COSE = zeros(size(powerSeg_ENV_COSIE,2),4);
bscEstimate_ENV_COSE(:,3) = powerSeg_ENV_COSIE(3,:);
bscEstimate_ENV_COSE(:,4) = powerSeg_ENV_COSIE(4,:);
bscEstimate_ENV_COSE(:,1) = (powerSeg_ENV_COSIE(1,:)./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;
bscEstimate_ENV_COSE(:,2) = ((powerSeg_ENV_COSIE(2,:))./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;

figure
hCohSpeckle = histogram(speckleCOSIE.thVector,'Normalization','probability');
binCentreCoh = hCohSpeckle.BinEdges(2:end)-hCohSpeckle.BinWidth;
pCoh = hCohSpeckle.Values;
close(gcf)

powerSeg_COH_WEIGHT = PDistWeighting(cohTest,binCentreCoh,pCoh,powf0,kWidth,oLap);
bscEstimate_COH_WEIGHT = zeros(size(powerSeg_COH_WEIGHT,2),4);
bscEstimate_COH_WEIGHT(:,3) = powerSeg_COH_WEIGHT(3,:);
bscEstimate_COH_WEIGHT(:,4) = powerSeg_COH_WEIGHT(4,:);
bscEstimate_COH_WEIGHT(:,1) = (powerSeg_COH_WEIGHT(1,:)./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;
bscEstimate_COH_WEIGHT(:,2) = ((powerSeg_COH_WEIGHT(2,:))./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;

figure
hEnvSpeckle = histogram(speckleSNRdata.snr,'Normalization','probability');
binCentreEnv = hEnvSpeckle.BinEdges(2:end)-hEnvSpeckle.BinWidth;
pEnv = hEnvSpeckle.Values;
close(gcf)


powerSeg_ENV_WEIGHT = PDistWeighting(snrEnv,binCentreEnv,pEnv,powf0,kWidth,oLap);
bscEstimate_ENV_WEIGHT = zeros(size(powerSeg_ENV_WEIGHT,2),4);
bscEstimate_ENV_WEIGHT(:,3) = powerSeg_ENV_WEIGHT(3,:);
bscEstimate_ENV_WEIGHT(:,4) = powerSeg_ENV_WEIGHT(4,:);
bscEstimate_ENV_WEIGHT(:,1) = (powerSeg_ENV_WEIGHT(1,:)./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;
bscEstimate_ENV_WEIGHT(:,2) = ((powerSeg_ENV_WEIGHT(2,:))./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;

%% Plot segmentation vs BSC 

bscEstimationSegFigure(bscSpeckleBf_test,bscSpeckleSTD_test,bscEstimate_COH_COSIE)
title('Coherence COSIE Segmentation')



%bmCohCOSIEparImage(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs,axIdxsCOH,powf0,segBool1_cluster,xBool,speckleCOSIE.pctSeg2(EMLidx),'Coherence')


bscEstimationSegFigure(bscSpeckleBf_test,bscSpeckleSTD_test,bscEstimate_ENV_COSE)
title('SNR COSIE Segmentation')

%bmCohCOSIEparImage(rayIdxs.*lWidth,yVals,bfImgData,depthIdx,axIdxs,axIdxsCOH,powf0,segBool2,xBool,speckleSNRdata.pctSeg2(EMLidx),'SNR')

%Pdist may have 0 segmentation percentage, so we'll plot with a seperate
%function to the COSIE results


bscEstimationWeightFigure(bscSpeckleBf_ref,bscSpeckleSTD_ref,bscEstimate_ENV_WEIGHT,bscEstimate_COH_WEIGHT)


save([saveDir_SEG,'\SegResults.mat'],'bscSpeckleBf_ref','bscSpeckleSTD_ref','bscEstimate_COH_COSIE',...
    'bscEstimate_ENV_COSE','bscEstimate_ENV_WEIGHT','bscEstimate_COH_WEIGHT')


