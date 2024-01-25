

% I will now make me and this fifty pence piece... disappear

clear

topDir = [uigetdir('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74_QuadInterp\','Select Analysis directory'),'\'];

sumIdx = 32;

fileNames = ls(topDir)

iImage = input('Image # : ');

for iFile = 3:size(fileNames,1)
    foldBool = isfolder([topDir,'\',fileNames(iFile,:)]);

    if foldBool
        [topDir,'\',fileNames(iFile,:)],
        ls([topDir,'\',fileNames(iFile,:)])
    end
end

load([topDir,'BFimgData',num2str(iImage),'.mat'])


adaptBool = 1;%input('Adaptive grid sizing? : ');

if adaptBool 
    adaptStr = '_adaptive';
else
    adaptStr = '';
end

speckleDir = ['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74_QuadInterp\Speckle\'];

vsxParams = load([topDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

depthSelect = 25;%input('Depth of interest (mm) : ');

speckleDir = [speckleDir,'Z',num2str(depthSelect),'\'];
topDir = [topDir,'Z',num2str(depthSelect),'\'];
load([topDir,'COSIEinput',num2str(iImage),'.mat'])


tgcBool = isequal(vsxParams.TGC,vsxParams2.TGC)

if ~tgcBool
    error('TGC altered between Speckle and test image')
end




%% 1. Set up variables 

samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample  + 1;

yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;

[~, depthIdx] = min(abs(yVals - depthSelect*1e-3));
kLength_BSC_samples = 120;
axIdxs = depthIdx-round(kLength_BSC_samples/2) : depthIdx + round(kLength_BSC_samples/2) -1 ;

fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);
 
powf0 = abs(spectAll(:,nF));

%sumIdx = size(cohAll,2);


%% 2. Compute segmentation based on EML

imgDir = [topDir,'\SegResults_',num2str(iImage),adaptStr,'\SumIdx_',num2str(sumIdx)];

if ~exist(imgDir)
    mkdir(imgDir)
end

speckleCOSIE = load([speckleDir,'\COSIEoutput',adaptStr,'\COSIEoutput',num2str(sumIdx),'.mat']);
speckleEnvelope = load([speckleDir,'\envData_COSIE.mat']);


cohSum = sum(cohAll(:,1:sumIdx),2);


%This step is for interest only, just want to see how the EML's differ
%[bscSurface,segSurface,EML,pctSeg1,redEML,pctSeg2] = COSIE_adaptiveGrid(cohSum,powf0);

transpData = zeros(size(fullIM));
colorData = nan.*zeros(size(fullIM));

kWidth = 5; 
oLap = 0.8;

tGCV = linspace(0.3,1,size(fullIM,2));
tGCM = repmat(tGCV,size(fullIM,1),1);

iSegPct = 1;
    
%1 = keep  , 0 = segment
segBool = cohSum > speckleCOSIE.EML(1,iSegPct) & cohSum < speckleCOSIE.EML(2,iSegPct);

outSeg = (1-length(find(segBool))/length(segBool))*100;

%% 3. Plot BSC/b-mode image
lWidth = (Trans.ElementPos(2,1)-Trans.ElementPos(1,1))*lambda;

ax0 = axes;
B80 = bmode(iq',80);
imagesc(ax0,lWidth.*(1:128),yVals,B80);
hold on 
plot(ax0,lWidth.*[1 size(fullIM,1)],yVals(depthIdx).*[1 1],'-','color','red')
plot(ax0,lWidth.*120.*[1 1], [yVals(axIdxs(1)) yVals(axIdxs(end))],'-','color','red')
xlabel('Lateral Position (m)')
axis equal
axis tight
ylabel('Axial Position (m)')
colormap(ax0,'gray') 
fName1 =[imgDir,'/Bmode'];
savefig(fName1)
saveas(gcf,fName1,'png')

ax1 = axes;
imagesc(ax1,lWidth.*(1:128),yVals,B80);
hold on 
plot(ax1,lWidth.*[1 size(fullIM,1)],yVals(depthIdx).*[1 1],'-','color','red')
plot(ax1,lWidth.*120.*[1 1], [yVals(axIdxs(1)) yVals(axIdxs(end))],'-','color','red')
xlabel('Lateral Position (m)')
axis equal
ylabel('Axial Position (m)')

ax2 = axes;
powfLong = repmat(powf0,1,size(axIdxs,2));

hideVectorLong = repmat(segBool,1,size(axIdxs,2));

colorData(rayIdxs,axIdxs) = log(powfLong); 

transpData(rayIdxs,axIdxs) = 0.6.*hideVectorLong; 

imagesc(ax2 ,lWidth.*(1:128) , yVals, colorData','AlphaData',transpData')    
cB = colorbar;
cB.Label.String = 'Log(\kappa)';
cB.Label.FontSize= 20;
cB.Position(3) = 0.02;

linkaxes([ax1,ax2])

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray')
%colormap(ax2,'winter')
%ylim([0.04 0.09])

title(ax1,['Seg = ',num2str(outSeg),'%'])
axis tight 
axis equal

fname2 = [imgDir,'/ParametricImage_COH'];
savefig(fname2)
saveas(gcf,fname2,'png')


cohValues2 = sum(RMat(:,:,1:sumIdx),3);
F = griddedInterpolant(xQcoh,yQcoh,cohValues2);
CoherenceQ = F(xQbm,yQbm);

figure 
imagesc(lWidth.*(1:128) , yVals, CoherenceQ');
xlabel('Lateral Position (m)')
ylabel('Axial Position (m)')
hold on 
plot(lWidth.*[1 size(fullIM,1)],yVals(depthIdx).*[1 1],'-','color','red')
plot(lWidth.*120.*[1 1], [yVals(axIdxs(1)) yVals(axIdxs(end))],'-','color','red')
colormap gray
axis tight 
axis equal

fname2 = [imgDir,'/CoherenceImage'];
savefig(fname2)
saveas(gcf,fname2,'png')

%% 5. Plot beam features at the focus

figure
subplot(2,2,1)
plot(lWidth.*rayIdxs,cohSum)
hold on
plot(lWidth.*[1 128],[speckleCOSIE.EML(2,1) speckleCOSIE.EML(2,1)],'-','Color','red')
plot(lWidth.*[1 128],[speckleCOSIE.EML(1,1) speckleCOSIE.EML(1,1)],'-','Color','red')
errorbar(lWidth.*120,mean(speckleCOSIE.thVector),std(speckleCOSIE.thVector),'.','Color','red')
ylim([-sumIdx/4 1.5*max(speckleCOSIE.EML(2,1))])
xlabel('Lateral Position')
ylabel('Summed Coherence')
xlim(lWidth.*[1 128])
title(['Coherence Seg.: ',num2str(outSeg)])


subplot(2,2,2)
hCohSpeckle = histogram(speckleCOSIE.thVector,'Normalization','probability');
binCentreCoh = hCohSpeckle.BinEdges(2:end)-hCohSpeckle.BinWidth;
pCoh = hCohSpeckle.Values;

pPoints = interp1(binCentreCoh,pCoh,cohSum);
pPoints(isnan(pPoints))= 0;
plot(lWidth.*rayIdxs,pPoints,'-','Color','k')
ylim([-0.01 0.15])
ylabel('P(Coherence)')
xlabel('Lateral Position')


subplot(2,2,3) 
focusEnv = abs(envelope(fullIM(rayIdxs,axIdxs)'));
mEnv = mean(focusEnv,1);
sEnv = std(focusEnv,0,1);
snrEnv = mEnv./sEnv;

segBool = snrEnv > speckleEnvelope.EML(1,iSegPct) & snrEnv < speckleEnvelope.EML(2,iSegPct);
outSeg = (1-length(find(segBool))/length(segBool))*100;


plot(lWidth.*rayIdxs,snrEnv,'-')
hold on 
plot(lWidth.*[1 128],[speckleEnvelope.EML(2,1) speckleEnvelope.EML(2,1)],'-','Color','red')
plot(lWidth.*[1 128],[speckleEnvelope.EML(1,1) speckleEnvelope.EML(1,1)],'-','Color','red')
plot(lWidth.*[1 size(fullIM,1)],1.91.*[1 1],'-','Color','red')
ylabel('SNR')
xlabel('Lateral Position')
xlim(lWidth.*[1 128])
title(['SNR Seg.: ',num2str(outSeg)])
errorbar(lWidth.*120,mean(snrEnv),std(snrEnv),'.','Color','red')

subplot(2,2,4)
hEnvSpeckle = histogram(speckleEnvelope.snr,'Normalization','probability');
binCentreEnv = hEnvSpeckle.BinEdges(2:end)-hEnvSpeckle.BinWidth;
pEnv = hEnvSpeckle.Values;

pPoints = interp1(binCentreEnv,pEnv,snrEnv);
pPoints(isnan(pPoints))= 0;
plot(lWidth.*rayIdxs,pPoints,'-','Color','k')
ylim([-0.01 0.15])
ylabel('P(SNR)')


fname2 = [imgDir,'\FocLines'];
savefig(fname2)
saveas(gcf,fname2,'png')




%% 5. Repeat 4. for more seg pctgs 
outSegPrev = 0; 

kWidth = 5;
oLap = 0.8;
kIdxs = idxClustering(1:size(spectAll,1),kWidth,oLap);
nKernels = size(kIdxs,2);
bscValuesIK = zeros(1,nKernels);
stdValuesIK = zeros(1,nKernels);

nPossibleKernels = size(idxClustering(1:length(cohSum),kWidth,oLap),2);


for iKernel = 1:nKernels
    bscValuesIK(iKernel) = abs(mean(spectAll(kIdxs{iKernel},nF)));
end


powf0 = spectAll(:,nF);

bscEstimate_COH = COVsegmentation(cohSum,speckleCOSIE.EML,powf0,kWidth,oLap,nPossibleKernels);
bscEstimate_ENV = COVsegmentation(snrEnv,speckleEnvelope.EML,powf0,kWidth,oLap,nPossibleKernels);


bscEstimate_weighted_COH = COVWeighting(cohSum,binCentreCoh,pCoh,powf0,kWidth,oLap);
bscEstimate_weighted_ENV = COVWeighting(snrEnv,binCentreEnv,pEnv,powf0,kWidth,oLap);


normDistData =  normrnd(zeros(1,1e5),ones(1,1e5));

figure
plot([std(speckleCOSIE.thVector)/mean(speckleCOSIE.thVector), skewness(speckleCOSIE.thVector) ,kurtosis(speckleCOSIE.thVector)],'o','MarkerFaceColor','blue')
hold on 
plot([std(speckleEnvelope.snr)/mean(speckleEnvelope.snr), skewness(speckleEnvelope.snr) ,kurtosis(speckleEnvelope.snr)],'o','MarkerFaceColor','red')
plot([nan, skewness(normDistData) ,kurtosis(normDistData)],'o','Color','white','MarkerFaceColor','k')
xticks([1,2,3])
xlim([0 4])
xticklabels({'C.O.V.','Skew','Kurtosis'})
ylabel('Value')
l = legend({'Coherence','Texture','Normal distn.'});
l.Position = [0.1482 0.7849 0.2185 0.1214];
sgtitle('Statistics of Coherence and Texture Distributions')


%% 6. Plots covariance/kurtosis data

%kurtosis/cov switch , 1 = kurtosis, 2 = COV
kCOVswitch = 2;

tString = 'BSC: Coherence analysis (COSIE)';
plotOptions(bscEstimate_COH,tString,kCOVswitch)
fname3 = [imgDir,'\VariabilityCOSIEcoh',];
savefig(fname3)
saveas(gcf,fname3,'png')


tString = 'BSC: SNR analysis (COSIE)';
plotOptions(bscEstimate_weighted_COH,tString,kCOVswitch)
fname4= [imgDir,'\VariabilityCOSIEsnr',];
savefig(fname4)
saveas(gcf,fname4,'png')


tString = 'BSC: Coherence analysis (weighting)';
plotOptions(bscEstimate_ENV,tString,kCOVswitch)
fname5 = [imgDir,'\VariabilityWeightcoh',];
savefig(fname5)
saveas(gcf,fname5,'png')


tString = 'BSC: SNR analysis (weighting)';
plotOptions(bscEstimate_weighted_ENV,tString,kCOVswitch)
fname6 = [imgDir,'\VariabilityWeightsnr',];
savefig(fname6)
saveas(gcf,fname6,'png')


%% 7. final boss figures

figure 
plot(bscEstimate_COH(3,:),bscEstimate_COH(2,:)./bscEstimate_COH(1,:),'-o','Color','k','MarkerFaceColor','k')
hold on 
plot(bscEstimate_weighted_COH(3,:),bscEstimate_weighted_COH(2,:)./bscEstimate_weighted_COH(1,:),'-sq','Color','k','MarkerFaceColor','k')
plot(bscEstimate_ENV(3,:),bscEstimate_ENV(2,:)./bscEstimate_ENV(1,:),'-o','Color','red','MarkerFaceColor','red')
plot(bscEstimate_weighted_ENV(3,:),bscEstimate_weighted_ENV(2,:)./bscEstimate_weighted_ENV(1,:),'-sq','Color','red','MarkerFaceColor','red')
xlim padded
xlabel('Seg %')
ylabel('C.O.V of BSC')
legend({'COSIE-COH','Weighted-COH','COSIE-ENV','Weighted-ENV'})
fname7 = [imgDir,'\COVfigure',];
savefig(fname7)
saveas(gcf,fname7,'png')


figure 
plot(bscEstimate_COH(3,:),bscEstimate_COH(5,:),'-o','Color','k','MarkerFaceColor','k')
hold on 
plot(bscEstimate_weighted_COH(3,:),bscEstimate_weighted_COH(5,:),'-sq','Color','k','MarkerFaceColor','k')
plot(bscEstimate_ENV(3,:),bscEstimate_ENV(5,:),'-o','Color','red','MarkerFaceColor','red')
plot(bscEstimate_weighted_ENV(3,:),bscEstimate_weighted_ENV(5,:),'-sq','Color','red','MarkerFaceColor','red')
xlim padded
xlabel('Seg %')
ylabel('Kurtosis of BSC')
legend({'COSIE-COH','Weighted-COH','COSIE-ENV','Weighted-ENV'})
fname8 = [imgDir,'\KurtosisFigure',];
savefig(fname8)
saveas(gcf,fname8,'png')


save([imgDir,'\StatAttack.mat'],'bscEstimate_COH','bscEstimate_weighted_COH','bscEstimate_ENV','bscEstimate_weighted_ENV')

