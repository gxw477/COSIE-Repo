


clear 
close all

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver_NHV_NTGC\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\G218L74_1\';

%load verasonics param's
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

iImage =  input('Which image ? : ');
1
bfImgData = load([testDir,'BFimgData',num2str(iImage),'.mat']);


kLength_BSC_samples = 120;

fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
fVals = (0:kLength_BSC_samples-1).*fs/(kLength_BSC_samples);
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);


%ground truth BSC values: reference
mu0_ref = 3.25e-4/(3^3.6);
bscSpeckleBf_ref = mu0_ref * vsxParams.Trans.frequency^3.6; 
bscSpeckleSTD_ref = 0.2*bscSpeckleBf_ref;

%ground truth BSC values: test
mu0_test = 3.86e-4/(3^3.5);
bscSpeckleBf_test = mu0_test * vsxParams.Trans.frequency^3.6; 
bscSpeckleSTD_test = 0.2*bscSpeckleBf_test;

%Calculate edge correction factor
[mEdgeSpectSpeckle , yEdgeSpeckle] = edgeDetectionMulti(speckleDir,1);
[mEdgeSpectTest ,yEdgeTest] = edgeDetectionMulti(testDir,0);
%[mEdgeSpectTest ,yEdgeTest] = edgeDetectionSingle(bfImgData.fullIM,bfImgData.fs,bfImgData.yVals,10,1540,bfImgData.Trans.frequency*1e6,17:111);

edgecorr = (mean(mEdgeSpectSpeckle)/mean(mEdgeSpectTest))^2;

%edgecorr =   46.6590;

bm = bmode(bfImgData.iq',80);  
[~ , edge] = (max(bm,[],1));
edgeYVal = bfImgData.yVals(round(mean(edge)));

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

EMLidx =1;
kWidth = 5; 
oLap = 0.8;

COSIEall = cell(1,length(depthVals));
powAll = zeros(length(depthVals), length(rayIdxs));
cohAll = zeros(length(depthVals), length(rayIdxs));
segBoolAll = zeros(length(depthVals), length(rayIdxs));

bscAll_COSIE = zeros(length(depthVals),2);
bscAll_UNSEG = zeros(length(depthVals),2);



for iDepth = 1:length(depthVals)

    depthSelect = depthVals(iDepth);

    dataDir = [testDir,wName,'\Z',num2str(depthSelect),'\'];
    % loop this over depths
    cInput = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);
    speckleDir2 = [speckleDir,wName,'\Z',num2str(depthSelect),'\'];  
    speckleCOSIE = load([speckleDir2,'\COSIEoutput',adaptStr,num2str(kLength_COH_samples),'\COSIEoutput',num2str(sumIdx),'.mat']);
    
    
    %centre frequency of test data
    powf0 = abs(cInput.spectAll(:,nF));

    attTest  = [0.579 , 0.955].*vsxParams.Trans.frequency;
  
    dzT = depthSelect*0.1 - edgeYVal*100;
    attTest_DB = 2*(dzT)*attTest(1);
    attComp_Test =   10^(attTest_DB/10);

    specklePOWER_MEAN = mean(speckleCOSIE.powf0);
   
    
    bscConvert = (1./specklePOWER_MEAN) * bscSpeckleBf_ref *edgecorr*attComp_Test;
    powAll(iDepth,:) = powf0.*bscConvert;
    
    RMatTest = zeros(length(rayIdxs),sumIdx);
    
    [~ , cohTestKernelIdx] = min(abs(depthSelect*1e-3 - bfImgData.yVals));
    axIdxsCOH = (cohTestKernelIdx- kLength_COH_samples/2):(cohTestKernelIdx+kLength_COH_samples/2-1);

    for iLine = 1:length(rayIdxs)
       RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(rayIdxs(iLine),axIdxsCOH,:)));
    end

    cohTest = sum(RMatTest(:,1:sumIdx),2);
    cohAll(iDepth,:) = cohTest;

    segBool = cohTest > speckleCOSIE.redEML(1,EMLidx) & cohTest < speckleCOSIE.redEML(2,EMLidx);
  
    segBool1_cluster = ismember(rayIdxs, unique(cell2mat(idxClustering(rayIdxs(segBool),kWidth,oLap))))';

    powerSeg_COH_COSIE = COVsegmentation_sK(cohTest,speckleCOSIE.redEML,powAll(iDepth,:),kWidth,oLap);


    bscAll_UNSEG(iDepth,:) = [powerSeg_COH_COSIE(1,1),powerSeg_COH_COSIE(2,1)];
    bscAll_COSIE(iDepth,:) = [powerSeg_COH_COSIE(1,EMLidx+1),powerSeg_COH_COSIE(2,EMLidx+1)];
    
    outSeg = 100*(1-length(find(segBool1_cluster))/length(segBool1_cluster));

    speckleScore(iDepth) = 100 - (outSeg-speckleCOSIE.pctSeg2(EMLidx));
    
    segPct(iDepth) = 100*(length(segBool)-length(find(segBool)))/length(segBool);

    segBoolAll(iDepth,:) = segBool1_cluster;

    COSIEall{iDepth} = speckleCOSIE;

end



samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample  + 1;
yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;

lambda = 1540/vsxParams.Trans.frequency*1e-6;
lWidth = lambda*(vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1)) ;
xVals = (1:128).*lWidth;


MultipleParametricCohImage(xVals.*1e3,yVals.*1e3,rayIdxs,bfImgData,depthVals,powAll,segBoolAll,kLength_BSC_samples,kLength_COH_samples)


bscSpeckleBf_ref = bscSpeckleBf_ref*1e2;
bscSpeckleSTD_ref = bscSpeckleSTD_ref*1e2;

bscSpeckleBf_test = bscSpeckleBf_test*1e2;
bscSpeckleSTD_test = bscSpeckleSTD_test*1e2;

%% Calculate depth1 features 


powAll2 = powAll;
powAll2(find(~segBoolAll)) = nan;

meanBSC_COSIE = mean(powAll2,2,'omitnan');
errBSC_COSIE = std(powAll2,[],2,'omitnan');


meanBSC_unseg = mean(powAll,2,'omitnan');
errBSC_unseg = std(powAll,[],2,'omitnan');

maxBSCall = max([meanBSC_COSIE(:)+errBSC_COSIE(:);meanBSC_unseg(:)+errBSC_unseg(:)]);

figure
errorbar(depthVals-0.5,meanBSC_COSIE,errBSC_COSIE,'.','MarkerSize',10)
hold on 
errorbar(depthVals+0.5,meanBSC_unseg,errBSC_unseg,'.','MarkerSize',10)
plot([min(depthVals)-5 max(depthVals)+5],bscSpeckleBf_test.*[1 1],'k-','MarkerSize',10)
errorbar(min(depthVals)-5,bscSpeckleBf_test, bscSpeckleSTD_ref,'k.','MarkerSize',10)
xticks(depthVals)
xlabel('Depth (mm)')
ylabel('BSC (m^{-1}sr^{-1})')
%title('Full Width averaging')
ax = gca;
ax.FontSize = 20;
xlim padded

ylim([0 1.2])
yticks(ax,[0:0.2:1.2])
%ylim padded

ax1 = gca;

figure
errorbar(depthVals-0.5,bscAll_COSIE(:,1),bscAll_COSIE(:,2),'.','MarkerSize',10)
hold on 
errorbar(depthVals+0.5,bscAll_UNSEG(:,1),bscAll_UNSEG(:,2),'.','MarkerSize',10)
plot([min(depthVals)-5 max(depthVals)+5],bscSpeckleBf_test.*[1 1],'k-','MarkerSize',10)
errorbar(min(depthVals)-5,bscSpeckleBf_test, bscSpeckleSTD_ref,'k.','MarkerSize',10)
xticks(depthVals)
set(gca,'FontSize', 20);

xlabel('Depth (mm)')
ylabel('BSC (cm^{-1}sr^{-1})')
title('Kernel averaging')

xlim(ax1.XLim)
ylim(ax1.YLim)


figure
plot(depthVals, speckleScore,'kx','MarkerSize',10)
xlim([7.55 47.45])
xticks(ax1.XTick)
xlabel('Depth (mm)')
ylabel('Speckle Score')
set(gca,'FontSize', 20);

xlim padded
ylim padded


errVals_COSIE = 100.*abs(bscAll_COSIE(:,1)- bscSpeckleBf_test)./bscSpeckleBf_test;
errVals_Unseg = 100.*abs(bscAll_UNSEG(:,1)- bscSpeckleBf_test)./bscSpeckleBf_test;

mError = max([errVals_COSIE(:);errVals_Unseg(:)]);

r25 = ceil(mError/25)*25;

cMatrix = colororder;
markerSizeMatrix = 8.*(100-segPct)/100;

figure
hold on 

for iPlotMarker = 1:length(markerSizeMatrix(:))

    plot(depthVals(iPlotMarker)-0.5, errVals_COSIE(iPlotMarker) ,'o','MarkerSize',markerSizeMatrix(iPlotMarker),'Color',cMatrix(1,:),'MarkerFaceColor',cMatrix(1,:))
    
end

plot(depthVals+0.5, errVals_Unseg,'o','MarkerSize', 6,'MarkerFaceColor',cMatrix(2,:))
xlabel('Depth (mm)')
xlim([7.55 47.45])
xticks(ax1.XTick)
ylabel('BSC error (%)')
%legend({'COSIE','Unsegmented'})
set(gca,'FontSize', 20);

xlim padded
ylim([0 r25])
yticks(0:25:r25)




