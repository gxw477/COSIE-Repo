

clear 
close all 

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QA\';

%load verasonics param's
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

sumIdx = 33;
iImage = input('Which Image ? :');
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
speckleCOSIE = load([speckleDir2,'\COSIEoutput',adaptStr,'_smallK\COSIEoutput',num2str(sumIdx),'.mat']);

dataDir = [testDir,wName,'\Z',num2str(depthSelect),'\'];

%contains spectAll 
powerInputStruct = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);


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
powf0 = abs(powerInputStruct.spectAll(:,nF));
bscKlength = length(powerInputStruct.spectAll(1,:));

%Calculate edge correction factor
%[mEdgeSpectSpeckle ,yEdgeSpeckle] = edgeDetectionMulti(speckleDir);
%[mEdgeSpectTest ,yEdgeTest] = edgeDetectionMulti(testDir);


edgecorr =   46.6590;%(mean(mEdgeSpectSpeckle)/mean(mEdgeSpectTest))


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

%how many coherence kernels fit in a BSC kernel? 
nCohKernpBSC = length(powerInputStruct.spectAll(1,:))/bfImgData.kLength;

if nCohKernpBSC/2 ~= round(nCohKernpBSC/2)
    error('Sort out kernel lengths')
end

[~ , cohTestKernelIdx] = min(abs(depthSelect*1e-3 - bfImgData.kY));
cohTest = sum(bfImgData.RMat(xBool,cohTestKernelIdx-nCohKernpBSC/2 : cohTestKernelIdx + nCohKernpBSC/2-1 ,1:sumIdx) ,3);



%% Generate parametric image using first segmentation point on EML

%The maximum percentage of rejected kernels we will accept for the BSC calc
rejectPCT = 90;
rejectThresh = nCohKernpBSC - 30/nCohKernpBSC;

segBool = cohTest > speckleCOSIE.EML(1,1) & cohTest < speckleCOSIE.EML(2,1);
segBool2 = sum(segBool,2);
segBool3 = segBool2 > rejectThresh;

bmCohCOSIEparImage((1:128).*lWidth,yVals,bfImgData,depthIdx,axIdxs,powf0,segBool3,xBool)

figure
plot(xVals(xBool),cohTest)
hold on
plot(lWidth.*[1 128],[speckleCOSIE.EML(2,1) speckleCOSIE.EML(2,1)],'-','Color','red')
plot(lWidth.*[1 128],[speckleCOSIE.EML(1,1) speckleCOSIE.EML(1,1)],'-','Color','red')
errorbar(lWidth.*120,mean(speckleCOSIE.thVector),std(speckleCOSIE.thVector),'.','Color','red')
ylim([-sumIdx/4 1.5*max(speckleCOSIE.EML(2,1))])
xlabel('Lateral Position')
ylabel('Summed Coherence')
xlim(lWidth.*[1 128])



%%

nEMLpoints = size(speckleCOSIE.EML,2);
%bscCOSIE = zeros(2,nEMLpoints+1);

kWidth = 5; 
oLap = 0.8;

powerSeg = COVsegmentation_mK(cohTest,speckleCOSIE.EML,(powf0),kWidth,oLap,rejectThresh);

bscEstimate = zeros(size(powerSeg,2),4);
bscEstimate(:,3) = powerSeg(3,:);
bscEstimate(:,4) = powerSeg(4,:);

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

[maxEstimate,maxEstimateIdx] = max(bscEstimate(1:end,1));
yLim2 = maxEstimate + bscEstimate(maxEstimateIdx,2)*1.1;

figure 
tL = tiledlayout(1,1);
ax1 = axes(tL);

errorbar(ax1,-10, bscSpeckleBf,bscSpeckleSTD,'r.')
hold on 
plot(ax1,[-10 bscEstimate(end,3)+10],bscSpeckleBf.*[1 1],'r-.')
errorbar(ax1,bscEstimate(1,3),bscEstimate(1,1),bscEstimate(1,2),'k.')
errorbar(ax1,bscEstimate(2:end,3),bscEstimate(2:end,1),bscEstimate(2:end,2),'k.')
xlabel('Segmentation (%)')
ylabel('BSC (m^{-1}sr^{-1})')
xticks(ax1,[0:10:80])
xlim(ax1,[-20 bscEstimate(end,3)+10])
ylim(ax1,[0 yLim2])

ax2 = axes(tL);
plot(ax2,bscEstimate(2:end,3),bscEstimate(2:end,1),'k.')
ylim(ax2,[0 yLim2])
xlim(ax2,[-20 bscEstimate(end,3)+10])
xticks(ax2,unique(bscEstimate(:,3)))
yticks(ax2,[])
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'left';
set(ax2,'YColor','k')
set(ax2,'XColor','k')

xticklabels(ax2,num2cell(round(bscEstimate(:,4))))
set(ax2,'FontSize',8)
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';


saveDir_SEG = [testDir,wName,'\Z',num2str(depthSelect),'\SegResults_',num2str(iImage),adaptStr,'_smallK1\SumIdx_',num2str(sumIdx),'\'];
  
if ~exist(saveDir_SEG,'dir')
    mkdir(saveDir_SEG)
end

%saveas(gcf,[saveDir_SEG,'ParametricImage_COH'])

%savefig(gcf,['BSC_EstimateFigure/S_',num2str(sumIdx)])
%saveas(gcf,['BSC_EstimateFigure/S_',num2str(sumIdx),'.jpg'])


%% more parametrics? 

contBool = input('More maps? ');

if contBool

    saveDirJPGS = [testDir,wName,'\Z',num2str(depthSelect),'\Jpgs\'];
    
    if ~exist(saveDirJPGS,'dir')
        mkdir(saveDirJPGS)
    end
    
    v = VideoWriter([saveDirJPGS,'myFile.mp4']);
    v.FrameRate = 1.5;
    open(v)
    
    for idxEML = 1:55

        segBool = cohTest > speckleCOSIE.EML(1,idxEML) & cohTest < speckleCOSIE.EML(2,idxEML);
        segBool2 = sum(segBool,2);
        segBool3 = segBool2 > rejectThresh;

        bmCohCOSIEparImage_EML_video((1:128).*lWidth,yVals,bfImgData,depthIdx,axIdxs,powf0,segBool3,xBool,speckleCOSIE,cohTest,sumIdx,idxEML)
        saveas(gcf,[saveDirJPGS,'idx_',num2str(idxEML),'.jpg']);
        A = imread([saveDirJPGS,'idx_',num2str(idxEML),'.jpg']);
        writeVideo(v,A);
        close all
    
    end
    
    close(v)

end


%% Calculations for segmentations vs BSC estimate with both methods and both input vbls

close all

speckleSNRdata= load([ speckleDir2 , 'envData_COSIE' ]) ;
qaSNRdata = load([testDir,wName,'\Z',num2str(depthSelect),'\EnvStats',num2str(iImage),'.mat']);

testEnvelope = load([dataDir,'\EnvStats',num2str(iImage),'.mat']);
    


powerSeg_COH_COSIE = COVsegmentation(cohTest,speckleCOSIE.EML,powf0,kWidth,oLap,rejectThresh);
bscEstimate_COH_COSIE = zeros(size(powerSeg_COH_COSIE,2),4);
bscEstimate_COH_COSIE(:,3) = powerSeg_COH_COSIE(3,:);
bscEstimate_COH_COSIE(:,4) = powerSeg_COH_COSIE(4,:);
bscEstimate_COH_COSIE(:,1) = (powerSeg_COH_COSIE(1,:)./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_COH_COSIE(:,2) = ((powerSeg_COH_COSIE(2,:))./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_COH_COSIE(:,5) = powerSeg_COH_COSIE(5,:);
bscEstimate_COH_COSIE(:,5) = powerSeg_COH_COSIE(6,:);


snrEnv = testEnvelope.envMean./testEnvelope.envStd;
powerSeg_ENV_COSIE = COVsegmentation(snrEnv,speckleSNRdata.EML,powf0,kWidth,oLap,rejectThresh);
bscEstimate_ENV_COSIE = zeros(size(powerSeg_ENV_COSIE,2),4);
bscEstimate_ENV_COSIE(:,3) = powerSeg_ENV_COSIE(3,:);
bscEstimate_ENV_COSIE(:,4) = powerSeg_ENV_COSIE(4,:);
bscEstimate_ENV_COSIE(:,1) = (powerSeg_ENV_COSIE(1,:)./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_ENV_COSIE(:,2) = ((powerSeg_ENV_COSIE(2,:))./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_ENV_COSIE(:,5) = powerSeg_ENV_COSIE(5,:);
bscEstimate_ENV_COSIE(:,6) = powerSeg_ENV_COSIE(6,:);

figure
hCohSpeckle = histogram(speckleCOSIE.thVector,'Normalization','probability');
binCentreCoh = hCohSpeckle.BinEdges(2:end)-hCohSpeckle.BinWidth;
pCoh = hCohSpeckle.Values;
close(gcf)

powerSeg_COH_WEIGHT = PDistWeighting(cohTest,binCentreCoh,pCoh,powf0,kWidth,oLap);
bscEstimate_COH_WEIGHT = zeros(size(powerSeg_COH_WEIGHT,2),4);
bscEstimate_COH_WEIGHT(:,3) = powerSeg_COH_WEIGHT(3,:);
bscEstimate_COH_WEIGHT(:,4) = powerSeg_COH_WEIGHT(4,:);
bscEstimate_COH_WEIGHT(:,1) = (powerSeg_COH_WEIGHT(1,:)./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_COH_WEIGHT(:,2) = ((powerSeg_COH_WEIGHT(2,:))./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_COH_WEIGHT(:,5) = powerSeg_COH_WEIGHT(5,:);
bscEstimate_COH_WEIGHT(:,6) = powerSeg_COH_WEIGHT(6,:);

figure
hEnvSpeckle = histogram(speckleSNRdata.snr,'Normalization','probability');
binCentreEnv = hEnvSpeckle.BinEdges(2:end)-hEnvSpeckle.BinWidth;
pEnv = hEnvSpeckle.Values;
close(gcf)


powerSeg_ENV_WEIGHT = PDistWeighting(snrEnv,binCentreEnv,pEnv,powf0,kWidth,oLap);
bscEstimate_ENV_WEIGHT = zeros(size(powerSeg_ENV_WEIGHT,2),4);
bscEstimate_ENV_WEIGHT(:,3) = powerSeg_ENV_WEIGHT(3,:);
bscEstimate_ENV_WEIGHT(:,4) = powerSeg_ENV_WEIGHT(4,:);
bscEstimate_ENV_WEIGHT(:,1) = (powerSeg_ENV_WEIGHT(1,:)./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_ENV_WEIGHT(:,2) = ((powerSeg_ENV_WEIGHT(2,:))./specklePOWER_MEAN) * bscSpeckleBf *edgecorr*attComp_Test;
bscEstimate_ENV_WEIGHT(:,5) = powerSeg_ENV_WEIGHT(5,:);
bscEstimate_ENV_WEIGHT(:,6) = powerSeg_ENV_WEIGHT(6,:);


%% Plot segmentation vs BSC 

bscEstimationSegFigure(bscSpeckleBf,bscSpeckleSTD,bscEstimate_COH_COSIE)
title('Coherence COSIE Segmentation')

bscEstimationSegFigure(bscSpeckleBf,bscSpeckleSTD,bscEstimate_ENV_COSIE)
title('SNR COSIE Segmentation')

%Pdist may have 0 segmentation percentage, so we'll plot with a seperate
%function to the COSIE results

bscEstimationWeightFigure(bscSpeckleBf,bscSpeckleSTD,bscEstimate_ENV_WEIGHT,bscEstimate_COH_WEIGHT)

%% plot statistics vs. Seg pct. 

statsBool = input('Statistics Plots ? : ');

if statsBool 

    plotOptions(bscEstimate_COH_COSIE','COSIE (coherence)',1)
    plotOptions(bscEstimate_COH_COSIE','COSIE (coherence)',2)
    
    plotOptions(bscEstimate_ENV_COSIE,'COSIE (env)',1)
    plotOptions(bscEstimate_ENV_COSIE,'COSIE (env)',2)
    
    plotOptions(bscEstimate_COH_WEIGHT,'Weighting (coherence)',1)
    plotOptions(bscEstimate_COH_WEIGHT,'Weighting (coherence)',2)
    
    plotOptions(bscEstimate_ENV_WEIGHT,'Weighting (env)',1)
    plotOptions(bscEstimate_ENV_WEIGHT,'Weighting (env)',2)

    normDistData =  normrnd(zeros(1,1e5),ones(1,1e5));
    
    figure
    plot([std(speckleCOSIE.thVector)/mean(speckleCOSIE.thVector), skewness(speckleCOSIE.thVector) ,kurtosis(speckleCOSIE.thVector)],'o','MarkerFaceColor','blue')
    hold on 
    plot([std(speckleSNRdata.snr)/mean(speckleSNRdata.snr), skewness(speckleSNRdata.snr) ,kurtosis(speckleSNRdata.snr)],'o','MarkerFaceColor','red')
    plot([nan, skewness(normDistData) ,kurtosis(normDistData)],'o','Color','white','MarkerFaceColor','k')
    xticks([1,2,3])
    xlim([0 4])
    xticklabels({'C.O.V.','Skew','Kurtosis'})
    ylabel('Value')
    l = legend({'Coherence','Texture','Normal distn.'});
    l.Position = [0.1482 0.7849 0.2185 0.1214];
    sgtitle('Statistics of Coherence and Texture Distributions')
  
    %savefig(fnameN)
    %saveas(gcf,fnameN,'png')
    
    

end

%% Save results 

save([saveDir_SEG,'\SegResults.mat'],'bscSpeckleBf','bscSpeckleSTD','bscEstimate_COH_COSIE',...
    'bscEstimate_ENV_COSIE','bscEstimate_ENV_WEIGHT','bscEstimate_COH_WEIGHT')









