

clear 
close all 

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1Dir\';

testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QA\';

sumIdx = 32;
iImage = 1;
bfImgData =  load([testDir,'BFimgData',num2str(iImage),'.mat']);


   
lWidth = (bfImgData.Trans.ElementPos(2,1)-bfImgData.Trans.ElementPos(1,1))*bfImgData.lambda;
   
depthSelect = 25;%input('Depth of interest (mm) : ');
    
adaptBool = 1;%input('Adaptive grid sizing? : ');

if adaptBool 
    adaptStr = '_adaptive';
else
    adaptStr = '';
end

vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

figure 
bm = bmode(bfImgData.iq',40);
imagesc(bm)

[~ , edge] = (max(bm,[],1));

edgeYVal =  bfImgData.yVals(round(mean(edge)));




wOption = 2;%input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');

if wOption == 1 
    wName = 'Rect';
elseif wOption == 2
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welch';
end

speckleDir2 = [speckleDir,wName,'\Z',num2str(depthSelect),'\'];
  
dataDir = [testDir,wName,'\Z',num2str(depthSelect),'\'];
load([dataDir,'COSIEinput',num2str(iImage),'.mat'])

samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample  + 1;

yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;

[~, depthIdx] = min(abs(yVals - depthSelect*1e-3));
kLength_BSC_samples = 120;
axIdxs = depthIdx-round(kLength_BSC_samples/2) : depthIdx + round(kLength_BSC_samples/2) -1 ;

fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);

powf0 = abs(spectAll(:,nF)).^2;


attTest  = [0.579 , 0.955].*vsxParams.Trans.frequency;
attSpeckle   = [0.524 , 0.09].*vsxParams.Trans.frequency;
attP     = [0.53  , 0.003].*vsxParams.Trans.frequency;


attTest_DB = 2*(depthSelect*0.1 - edgeYVal*100)*attTest(1);
attComp_Test =   10^(attTest_DB/10);

attSpeckle_DB = 2*(depthSelect*0.1 - edgeYVal*100)*attTest(1);
attComp_Speckle= 10^(attSpeckle_DB/10);

attSpeckle_DB2 = 2*(depthSelect*0.1 - edgeYVal*100)*attP(1);
attComp_Speckle2= 10^(attSpeckle_DB/10);

speckleCOSIE = load([speckleDir2,'\COSIEoutput',adaptStr,'\COSIEoutput',num2str(sumIdx),'.mat']);

xBool = 17:111;
kLengthPER = 10;

load('ElastPhtL74_1607/Img1Dir/edgeCorr.mat')
%edgecorr = (mean(mEdgeSpectSpeckle)/mEdgeSpect1)^2;


specklePOWER_MEAN = mean(speckleCOSIE.powf0.^2*attComp_Speckle);
speckle_STD  = std(speckleCOSIE.powf0.^2*attComp_Speckle);

testPOWER_MEAN = mean(powf0*attComp_Test);
testPOWER_STD = std(powf0*attComp_Test);

   

mu0 = 3.25e-4/(3^3.6);

bscSpeckleBf = mu0 * vsxParams.Trans.frequency^3.6; 
bscSpeckleSTD = 0.2*bscSpeckleBf;

bsc1 = (testPOWER_MEAN./specklePOWER_MEAN) * bscSpeckleBf;
bscErr1 = std((testPOWER_MEAN./specklePOWER_MEAN).^1 * bscSpeckleBf);
bscMean1 = mean(bsc1)*edgecorr;

figure 
errorbar(0, bscSpeckleBf,bscSpeckleSTD,'r.')
hold on 
errorbar(1, bscMean1,bscSpeckleSTD,'k.')


[~ , cohTestKernelIdx] = min(abs(depthSelect*1e-3 - bfImgData.kY));
cohTest = sum(bfImgData.RMat(xBool,cohTestKernelIdx,:) ,3);

nEMLpoints = size(speckleCOSIE.EML,2);

%bscCOSIE = zeros(2,nEMLpoints+1);

kWidth = 5; 
oLap = 0.8;

nPossibleKernels = size(idxClustering(1:length(cohTest),kWidth,oLap),2);    
bscEstimate      = COVsegmentation(cohTest,speckleCOSIE.EML,powf0,kWidth,oLap,nPossibleKernels)
    

figure 
errorbar(bscEstimate(3,:),bscEstimate(1,:),bscEstimate(2,:),'k-.')



    

