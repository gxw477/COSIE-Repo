

clear 
close all 

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QA\';

%load verasonics param's
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

sumIdx = 32;
iImage = 1;
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



speckleDir2 = [speckleDir,wName,'\Z',num2str(depthSelect),'\'];  
%load COSIE data
speckleCOSIE = load([speckleDir2,'\COSIEoutput',adaptStr,'\COSIEoutput',num2str(sumIdx),'.mat']);

dataDir = [testDir,wName,'\Z',num2str(depthSelect),'\'];
load([dataDir,'COSIEinput',num2str(iImage),'.mat'])


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
edgeYVal = bfImgData.yVals(round(mean(edge)));

%spectral calculations
fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);

%centre frequency of test data
powf0 = abs(spectAll(:,nF));

%Calculate edge correction factor
[mEdgeSpectSpeckle ,yEdgeSpeckle] = edgeDetectionMulti(speckleDir);
[mEdgeSpectTest ,yEdgeTest] = edgeDetectionMulti(testDir);
edgecorr = (mean(mEdgeSpectSpeckle)/mean(mEdgeSpectTest));


attTest  = [0.579 , 0.955].*vsxParams.Trans.frequency;
attSpeckle   = [0.524 , 0.09].*vsxParams.Trans.frequency;
attP     = [0.53  , 0.003].*vsxParams.Trans.frequency;

dzT = depthSelect*0.1 - edgeYVal*100;
attTest_DB = 2*(dzT)*attTest(1);
attComp_Test =   10^(attTest_DB/10);

%mean speckle power
specklePOWER_MEAN = mean(speckleCOSIE.powf0);
speckle_STD  = std(speckleCOSIE.powf0);

%mean test power
testPOWER_MEAN = mean(powf0*attComp_Test);
testPOWER_STD = std(powf0*attComp_Test);

%ground truth BSC calcs
mu0 = 3.25e-4/(3^3.6);
bscSpeckleBf = mu0 * vsxParams.Trans.frequency^3.6; 
bscSpeckleSTD = 0.2*bscSpeckleBf;

%regions of image with full aperture in effect
xBool = 17:111;

%calculate coherence properties
[~ , cohTestKernelIdx] = min(abs(depthSelect*1e-3 - bfImgData.kY));
cohTest = sum(bfImgData.RMat(xBool,cohTestKernelIdx,:) ,3);


%% Generate parametric image using first segmentation point on EML

segBool = cohTest > speckleCOSIE.EML(1,1) & cohTest < speckleCOSIE.EML(2,1);
bmCohCOSIEparImage((1:128).*lWidth,yVals,bfImgData,depthIdx,axIdxs,powf0,segBool,xBool)


%%

nEMLpoints = size(speckleCOSIE.EML,2);
%bscCOSIE = zeros(2,nEMLpoints+1);

kWidth = 5; 
oLap = 0.8;

nPossibleKernels = size(idxClustering(1:length(cohTest),kWidth,oLap),2);    
powerSeg      = COVsegmentation(cohTest,speckleCOSIE.EML,(powf0),kWidth,oLap,nPossibleKernels);

bscEstimate = zeros(size(powerSeg,2),3);
bscEstimate(:,3) = powerSeg(3,:);
bscEstimate(:,1) = (powerSeg(1,:)./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate(:,2) = ((powerSeg(2,:))./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;

errs = zeros(size(bscEstimate,1),6);
%[perrSi,perrSp,perralphaT,perralphaR,perrbscR , TOTAL] x nSegmentations

for i = 1:size(bscEstimate,1)

    errs(i,:) = errorPropagationBSC(powerSeg(1,i),powerSeg(2,i),specklePOWER_MEAN,speckle_STD,dzT,attTest(2),2e-2,attSpeckle(2),bscSpeckleBf,0.2*bscSpeckleBf,bscEstimate(i,1));
    %bscEstimate(i,3) = errs(i,6);

end

figure 
errorbar(-10, bscSpeckleBf,bscSpeckleSTD,'r.')
hold on 
errorbar(bscEstimate(1,3),bscEstimate(1,1), errs(1,6),'k.')
errorbar(bscEstimate(:,3),bscEstimate(:,1), errs(:,6),'k-.')


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


figure 
errorbar(-10, bscSpeckleBf,bscSpeckleSTD,'r.')
hold on 
errorbar(bscEstimate(1,3),bscEstimate(1,1),bscEstimate(1,2),'k.')
errorbar(bscEstimate(2:end,3),bscEstimate(2:end,1),bscEstimate(2:end,2),'k.')
xlabel('Segmentation (%)')
ylabel('BSC (m^{-1}sr^{-1})')
xticks([0:10:80])
xlim([-20 80])
%xlim([-1 2])



    

