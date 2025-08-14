
clear 
close all 

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';

%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_HV_HTGC\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_NHV_NTGC\';

planeDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Pref\';

%load verasonics param's2
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);


path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE\COSIE-Repo\Functions\')
path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE\COSIE-Repo\Scripts\')
path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\AttenuationTesting\AttenuationAlgorithm\AttenuationGUI')

iImage = 1;%input('Which Image ? : ');
%attStruct = load(['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_NHV_NTGC\AttDataInterp2\Frame',num2str(iImage),'.mat']);
attStruct = load([ 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\IDFFiltAverage.mat']);

sumIdx = 33;%input('Sum Idx: ');




imgDir = ['C:\Users\gwest\Documents\COSIE paper 1\EmmaLiver\Img',num2str(iImage),'\AttSeg\'];


%load test data
bfImgData = load([testDir,'BFimgData',num2str(iImage),'.mat']);
depthSelect = 25;%input('Depth of interest (mm) : ');
adaptBool = 1;%input('Adaptive grid sizing? : ');

if adaptBool 
    adaptStr = '_adaptive';
else
    adaptStr = '';
end

wOption = 4;%input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch \n 4 for Hanning: \n ');

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
cohKlength_Length = 0.5*cohKlength*vsxParams.lambda/vsxParams.sPerWaveInit;


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


%% Ground truth values

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
 
%input('Check which attenuation value you are using : ')
%load('C:\Users\gwest\Documents\COSIE paper 1\EmmaLiver\RejIdxs.mat')

skinThick = 2e-3;
fatThick = 2.5e-3; 
muscleThick = 9e-3;

%Emma Liver 
skin.att = [21.158	1]; %Np/m/MHz
skin.rho = [1109]; %kg/m3
skin.c = [1624]; %m/s

muscle.att = [7.166	0.600441504]; %Np/m/MHz
muscle.rho = [1090];
muscle.c = [1588];

fat.att = [9.3191	1] ;%Np/m/MHz
fat.rho = [911];
fat.c = [1477];

liver.att = [35.958];
liver.rho = [1079];


r1 = abs((skin.rho*skin.c - fat.rho*fat.c)/(skin.rho*skin.c + fat.rho*fat.c));
r2 = abs((muscle.rho*muscle.c - fat.rho*fat.c)/(muscle.rho*muscle.c + fat.rho*fat.c));

Np2dB = 20/log(10);

sAtt = skin.att(1)*Np2dB* vsxParams.Trans.frequency; %dB/m
mAtt = muscle.att(1)*Np2dB * vsxParams.Trans.frequency; %dB/m
fAtt = fat.att(1)*Np2dB* vsxParams.Trans.frequency; %dB/m

liverAtt = liver.att(1)* Np2dB; %dB/m 


attSubCut = 2*(skinThick*sAtt + fatThick*fAtt + muscleThick*mAtt); %dB

attSpeckle  = [0.524 , 0.9].*vsxParams.Trans.frequency;
%attP        = [0.53  , 0.003].*vsxParams.Trans.frequency;


%ground truth BSC values: reference 
mu0_ref = 3.86e-4/(3^3.5); %(cm^{-1}sr^{-1}) 
bscSpeckleBf_ref = mu0_ref * 100 * vsxParams.Trans.frequency^3.5 ; %(m^{-1}sr^{-1})
bscSpeckleSTD_ref = [0.8 1.2].*mu0_ref;

mubsTH = 0.1;
mubsTH_STD = [0.04 2.2];

%% Segmentation setup

%regions of image with full aperture in effect
xBool = 17:112;

rayIdxs = (1:128);    
rayIdxs2 = rayIdxs(xBool);  
EMLidx = 10;
kWidth = 5; 
oLap = 0.8;

allDepths = 15:5:50;

segBoolBIG1 = true(length(rayIdxs2),length(allDepths));%zeros(length(rayIdxs2),length(allDepths));
segBoolBIG2 = segBoolBIG1;

analFolder = [testDir,'\AttDataInterp2\'];
frameNames = ls([analFolder,'\Frame*']);
BAEdata = load([analFolder,frameNames(iImage,:)]);

nLines = length(rayIdxs2);
nKernelsAtt = size(attStruct.Filt.F,2);
nKernelsBSC = length(allDepths);

segBoolAtt1 = true(nLines,nKernelsAtt);
%segBoolAtt2 = zeros(size(segBoolAtt1));
%segBoolFull2 = true(length(rayIdxs2),length(bfImgData.yVals));

powf0_BIG =  zeros(nLines,nKernelsBSC);

axIdxs_BIG = [];

allThValsCOH = zeros(2,nKernelsBSC);
allThValsSNR = zeros(2,nKernelsBSC);

allCOH = zeros(nLines,nKernelsBSC);
allSNR = allCOH;
       
   
for iDepth = 1:length(allDepths)

    speckleDir2 = [speckleDir,wName,'\Z',num2str(allDepths(iDepth)),'\'];  
    dataDir = [testDir,wName,'\Z',num2str(allDepths(iDepth)),'\'];
    
    dzT = allDepths(iDepth)*0.1 - edgeYVal*100; %cm 

    [~ , cohTestKernelIdx] = min(abs(allDepths(iDepth)*1e-3 - bfImgData.yVals));

    RMatTest = zeros(length(xBool),sumIdx);
    axIdxs = (cohTestKernelIdx- cohKlength/2):(cohTestKernelIdx+cohKlength/2-1);
    
    for iLine = 1:length(xBool)
       RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(xBool(iLine),axIdxs,1:sumIdx)));
    end
    
    axIdxs_BIG = [axIdxs_BIG,axIdxs];

    %% Coh+power analysis

   
    %cInput = load([dataDir,'/Sum',num2str(sumIdx),'/COSIEinput',num2str(iImage),'.mat']);

    speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(cohKlength),'\COSIEoutput',num2str(sumIdx),'.mat']);
    cInput_Depth = load([dataDir,'\Sum',num2str(sumIdx),'\COSIEinput',num2str(iImage),'.mat']);

    powf0 = abs(cInput_Depth.spectAll(:,nF));
    specklePOWER_MEAN = mean(speckleCOSIE.powf0);

    convFactor = 1./specklePOWER_MEAN * bscSpeckleBf_ref *edgecorr;%*attComp_Test;

    bscEstimate = (powf0.*convFactor);
    powf0_BIG(:,iDepth) = bscEstimate;
    
    cohTest = sum(RMatTest(:,1:sumIdx),2);

    %Store coherence value and EML 
    %allThValsCOH(:,iDepth) = [speckleCOSIE.redEML(1,EMLidx),speckleCOSIE.redEML(2,EMLidx)];
    allCOH(:,iDepth) = cohTest;
    
    

    %% SNR analysis 

    qaSNRdata = load([testDir,wName,'\Z',num2str(allDepths(iDepth)),'\Sum',num2str(sumIdx),'\EnvStats',num2str(iImage),'.mat']);
    snrTest = qaSNRdata.envMean./qaSNRdata.envStd;
    
    %load COSIE data
    speckleSNRdata= load([speckleDir2 ,'envData_COSIE' ]) ;
  
    %allThValsSNR(:,iDepth) = [speckleSNRdata.redEML(1,EMLidx),speckleSNRdata.redEML(2,EMLidx)];
    allSNR(:,iDepth) = snrTest;


end


%% Attenuation correction calculation

%liverStart = fatThick + muscleThick + skinThick;
liverStart = 0.014;

analFolder = [testDir,'\AttDataInterp\'];
frameNames = ls([analFolder,'\Frame*']);
BAEdata = load([analFolder,frameNames(iImage,:)]);
%nLines = length(rayIdxs2);
%nKernels = size(IDF.Filt.F,2);

nEMLidx = 10;
idxFactor = 1; 
nDepth = length(allDepths);
allAtt_COH = zeros(nEMLidx,nDepth,2);
allAtt_SNR = zeros(nEMLidx,nDepth,2);
allAtt_unseg = zeros(nDepth,2);


for iEMLidx = 1:nEMLidx
    
    segBoolAttCOH = false(nLines,nDepth);
    segBoolAttSNR = segBoolAttCOH;

    %segBoolAttSNR = true(length(BAEdata.BAE.IDFres.x),1);
    %segBoolAttSNR = reshape(segBoolAttSNR,[nLines,nKernels]);
    

    for iDepth = 1:nDepth
    
        speckleDir2 = [speckleDir,wName,'\Z',num2str(allDepths(iDepth)),'\'];  
        speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(cohKlength),'\COSIEoutput',num2str(sumIdx),'.mat']);
   
        [~ , cohTestKernelIdx] = min(abs(allDepths(iDepth)*1e-3 - bfImgData.yVals));
        axIdxs = (cohTestKernelIdx- cohKlength/2):(cohTestKernelIdx+cohKlength/2-1);
       
        liverThick = allDepths(iDepth)*1e-3 - liverStart;
        
        
        %reproduce the coherence segmentation bool across 'n' depths (kernels
        %set within attenuation calculation) 
        segBoolCOH = allCOH(:,iDepth) > speckleCOSIE.redEML(1,(iEMLidx-1)*idxFactor+1) & cohTest < speckleCOSIE.redEML(2,(iEMLidx-1)*idxFactor+1);    
        segBoolCluster_COH = ismember(rayIdxs2, unique(cell2mat(idxClustering(rayIdxs2(segBoolCOH),kWidth,oLap))))';
        
        %which depths does the coherence segmentation apply to? 
        segBoolAttCOH(:,iDepth) = segBoolCluster_COH;

        if allDepths(iDepth) > 1e3*(liverStart + 5e-3) 
            
            attLiverStruct = attenuationAnalyse(bfImgData,liverStart*100,allDepths(iDepth)*0.1,segBoolAttCOH,attStruct,rayIdxs2,allDepths);
            close all
 
            attLiver = attLiverStruct.alpha*1e6*vsxParams.Trans.frequency;
            attLiverFilt = attLiverStruct.alpha_filt*1e6*vsxParams.Trans.frequency;
                
            attInLiver = 2*liverThick * attLiver
            attInLiverFilt = 2*liverThick * attLiverFilt 
            iDepth;
            
            corr_COH = attLiverStruct.corr;
            corrFilt_COH = attLiverStruct.corr_filt;
        else
            attInLiver = nan; 
            attInLiverFilt = nan;
            corr_COH = nan; 
            corrFilt_COH = nan;
        end

        speckleDir2 = [speckleDir,wName,'\Z',num2str(allDepths(iDepth)),'\'];  
        speckleCOSIE_SNR =  load([speckleDir2,'\envData_COSIE.mat']);
    
        %reproduce the coherence segmentation bool across 'n' depths (kernels
        %set within attenuation calculation) 
        segBoolSNR = allSNR(:,iDepth) > speckleCOSIE_SNR.redEML(1,(iEMLidx-1)*idxFactor+1) & allSNR(:,iDepth) < speckleCOSIE_SNR.redEML(2,(iEMLidx-1)*idxFactor+1);    
        segBoolCluster_SNR = ismember(rayIdxs2, unique(cell2mat(idxClustering(rayIdxs2(segBoolSNR),kWidth,oLap))))';
        
        %which depths does the coherence segmentation apply to? 
        segBoolAttSNR(:,iDepth) = segBoolCluster_SNR;
    
        
        if allDepths(iDepth) > 1e3*(liverStart + 5e-3) 
            
            attLiverStruct = attenuationAnalyse(bfImgData,liverStart*100,allDepths(iDepth)*0.1,segBoolAttSNR,attStruct,rayIdxs2,allDepths);
            close all
 
            attInLiver_SNR = attLiverStruct.alpha_filt*1e6*vsxParams.Trans.frequency;
                
            attInLiverFilt_SNR = 2*liverThick * attInLiver_SNR; 
            iDepth;
            
            corr_SNR = attLiverStruct.corr;
            corrFilt_SNR = attLiverStruct.corr_filt;
        else
            attInLiver_SNR = nan; 
            attInLiverFilt_SNR = nan;
            corr_SNR = nan; 
            corrFilt_SNR = nan;
        end
    
    
        tempNew(iEMLidx,iDepth,:) = [attInLiver, attInLiverFilt,2* liverThick * liverAtt];

        attTest_DB = attSubCut + attInLiver;
        attTest_DB_Filt= attSubCut + attInLiverFilt;
        
        attComp_Test =  10^(attTest_DB/10);
        attComp_Test_Filt =  10^(attTest_DB_Filt/10);
        
        allAtt_unseg(iDepth,:) = [attComp_Test,corr_COH];
        allAtt_COH(iEMLidx,iDepth,:) = [attComp_Test_Filt,corrFilt_COH];
    
        segBoolBIG_COH(iEMLidx,iDepth,:) = segBoolCluster_COH;


        attTest_DB_SNR = attSubCut + attInLiver_SNR;
        attTest_DB_Filt_SNR = attSubCut + attInLiverFilt_SNR;
        
        attComp_Test =  10^(attTest_DB_SNR/10);
        attComp_Test_Filt =  10^(attTest_DB_Filt_SNR/10);
        
        allAtt_SNR(iEMLidx,iDepth,:) = [attComp_Test_Filt,corrFilt_SNR];
    
        segBoolBIG_SNR(iEMLidx,iDepth,:) = segBoolCluster_SNR;

    end
end

%% Attenuation plotting


atttenuationPlotter(allDepths,allAtt_unseg,allAtt_COH,liverStart,liverAtt,attSubCut)

atttenuationPlotter(allDepths,allAtt_unseg,allAtt_SNR,liverStart,liverAtt,attSubCut)




close all



%% BSC values for all Depths, all EMLidxs 

powf0_BIGGER = repmat(powf0_BIG,[1,1,nEMLidx]);

powf0BIG_unseg = powf0_BIG.*allAtt_unseg(:,1)';


powf0BIG_COH = zeros(nEMLidx,nDepth,length(rayIdxs2));
powf0BIG_SNR = zeros(nEMLidx,nDepth,length(rayIdxs2));

close all

plot(allAtt_unseg(:,1))

for iEMLidx = 1:nEMLidx 
    for iDepth = 1:nDepth
        powf0BIG_COH(iEMLidx,iDepth,:) = powf0_BIG(:,iDepth).*ones(length(rayIdxs2),1)*allAtt_COH(iEMLidx,iDepth,1);
        powf0BIG_SNR(iEMLidx,iDepth,:) = powf0_BIG(:,iDepth).*ones(length(rayIdxs2),1)*allAtt_SNR(iEMLidx,iDepth,1);
   
    end

    subplot(1,2,1)
    plot(allDepths,allAtt_COH(iEMLidx,:,1),'-o')
    hold on 
    plot(allDepths,allAtt_SNR(iEMLidx,:,1),'-o')
    
    subplot(1,2,2)
    plot(allDepths,mean(powf0BIG_COH(iEMLidx,:,:),3),'-o')
    hold on 
    plot(allDepths,mean(powf0BIG_SNR(iEMLidx,:,:),3),'-o')
    

end


segPow1 = nan([size(powf0BIG_COH)]);
segPow2 = nan([size(powf0BIG_COH)]);
segPow3 = nan([size(powf0BIG_COH)]);
segPow4= nan([size(powf0BIG_COH)]);

segPow1(segBoolBIG_COH) = powf0BIG_COH(segBoolBIG_COH); 
segPow2(segBoolBIG_SNR) = powf0BIG_SNR(segBoolBIG_SNR); 

bsc1 = zeros(nEMLidx,nDepth,2);
bsc1(:,:,1) = mean(segPow1,3,'omitnan');
bsc1(:,:,2) = std(segPow1,1,3,'omitnan');


bsc2 = zeros(nEMLidx,nDepth,2);
bsc2(:,:,1) = mean(segPow2,3,'omitnan');
bsc2(:,:,2) = std(segPow2,1,3,'omitnan');

jitter = linspace(-2,2,nEMLidx);

greenValues = round(linspace(255,0,nEMLidx));
cmap = [255.*ones(nEMLidx,1),greenValues',zeros(nEMLidx,1)]./255;

mubsTH = 0.1;
mubsTH_STD = [0.04 2.2];
neg = mubsTH_STD(1);
pos = mubsTH_STD(2);

fillColor = [140, 222, 162]./256;

figure 
xShade = [-15 -15 110 110 ];
yShade = [ neg pos pos neg ];
area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',neg,'ShowBaseLine','off');
hold on 
plot([-15 110],mubsTH.*[1,1],'k-.','LineWidth',2)
plot([-15 110],pos.*[1 1],'-.','LineWidth',2,'Color',fillColor)
plot([-15 110],neg.*[1 1],'-.','LineWidth',2,'Color', fillColor)
plot([-15 -15],[neg pos],'-.','LineWidth',2,'Color',fillColor)
plot([110 110],[neg pos],'-.','LineWidth',2,'Color',fillColor)
set(gca,'FontSize',20)


for iEMLidx = 1:nEMLidx
    
    errorbar(jitter(iEMLidx) + allDepths,bsc1(iEMLidx,:,1),bsc1(iEMLidx,:,2),'.','Color',cmap(iEMLidx,:),'MarkerSize',10)

end
errorbar(allDepths,mean(powf0BIG_unseg,1),std(powf0BIG_unseg,1,1),'k.','MarkerSize',10)

%title('Coherence Segmentation')
xlim([15 55])
set(gca,'YScale','log')
yticks(10.^[-7:2:5])
ylim([1e-7 1e5])
xlabel('Depth (mm)')
ylabel('BSC (m^{-1}rad^{-1})')
saveas(gcf,[imgDir,'\BSCImage1_COH.fig'])
exportgraphics(gcf,[imgDir,'\BSCImage1_COH.pdf'])



figure
xShade = [-15 -15 110 110 ];
yShade = [ neg pos pos neg ];
area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',neg,'ShowBaseLine','off');
hold on 
plot([-15 110],mubsTH.*[1,1],'k-.','LineWidth',2)
plot([-15 110],pos.*[1 1],'-.','LineWidth',2,'Color',fillColor)
plot([-15 110],neg.*[1 1],'-.','LineWidth',2,'Color', fillColor)
plot([-15 -15],[neg pos],'-.','LineWidth',2,'Color',fillColor)
plot([110 110],[neg pos],'-.','LineWidth',2,'Color',fillColor)

set(gca,'FontSize',20)

for iEMLidx = 1:nEMLidx
    
    errorbar(jitter(iEMLidx) + allDepths,bsc2(iEMLidx,:,1),bsc2(iEMLidx,:,2),'.','Color',cmap(iEMLidx,:),'MarkerSize',10)
    hold on 
    
end

errorbar(allDepths,mean(powf0BIG_unseg,1),std(powf0BIG_unseg,1,1),'k.','MarkerSize',10)

%title('Coherence Segmentation')
xlim([15 55])
set(gca,'YScale','log')
yticks(10.^[-7:2:5])
ylim([1e-7 1e5])
xlabel('Depth (mm)')
ylabel('BSC (m^{-1}rad^{-1})')

saveas(gcf,[imgDir,'\BSCImage1_SNR.fig'])
exportgraphics(gcf,[imgDir,'\BSCImage1_SNR.pdf'])



%bsc2 = [mean(segPow2,,'omitnan')', std(segPow2,1,'omitnan')'];
%bsc3 = [mean(segPow3,1,'omitnan')', std(segPow3,1,'omitnan')'];
%bsc4 = [mean(segPow4,1,'omitnan')', std(segPow4,1,'omitnan')'];



%% Parametric images


[ax1,ax2,cB] = bmCohCOSIEparImage_mDepth(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs_BIG,kLength_BSC_samples,powf0_BIGCorr2,segBoolBIG1,xBool,speckleCOSIE.pctSeg2(EMLidx),'Coherence');
xlim([-15 15])
ax1.XTick = [-15:5:15];
ax1.YTick = [2,5:10:55];
ax1.YLim = [2 55];
set(ax1,'FontSize',14)
set(ax1,'FontWeight','normal')
%saveas(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage)])
%saveas(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage),'.jpg'])
%savefigPDF_Crop(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage)])


[ax3,ax4,cB2] = bmCohCOSIEparImage_mDepth(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs_BIG,kLength_BSC_samples,powf0_BIGCorr2,segBoolBIG2,xBool,speckleSNRdata.pctSeg2(EMLidx),'Coherence');
xlim([-15 15])
ax3.XTick = [-15 :5:15];
ax3.YTick = [2 ,5:10:55];
ax3.YLim = [2 55];
set(ax3,'FontSize',18)
set(ax3,'FontWeight','normal')
%saveas(gcf,[saveDir_SEG,'ParametricImage_SNR',num2str(iImage)])
%saveas(gcf,[saveDir_SEG,'ParametricImage_SNR',num2str(iImage),'.jpg'])
%savefigPDF_Crop(gcf,[saveDir_SEG,'ParametricImage_SNR',num2str(iImage)])





