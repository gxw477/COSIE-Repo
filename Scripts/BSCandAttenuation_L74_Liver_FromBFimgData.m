
clear 
close all 

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';

%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_HV_HTGC\';
%testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_NHV_NTGC\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiverReducedSet\';


planeDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Pref\';

%load verasonics params2
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);


path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE\COSIE-Repo\Functions\')
path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE\COSIE-Repo\Scripts\')
path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\AttenuationTesting\AttenuationAlgorithm\AttenuationGUI')

iImage = input('Which Image ? : ');
%attStruct = load(['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_NHV_NTGC\AttDataInterp2\Frame',num2str(iImage),'.mat']);
attStruct = load([ 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\IDFFiltAverage.mat']);

sumIdx = 33;%input('Sum Idx: ');
kWidth = 5; 
oLap = 0.8;


imgDir = ['C:\Users\gwest\Documents\COSIE paper 1\EmmaLiver\Img',num2str(iImage),'\AttSeg\'];

if ~exist(imgDir,'dir')
    mkdir(imgDir)
end


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
%EMLidx = 5;
kWidth = 5; 
oLap = 0.8;

allDepths = 15:5:50;

segBoolBIG1 = true(length(rayIdxs2),length(allDepths));%zeros(length(rayIdxs2),length(allDepths));
segBoolBIG2 = segBoolBIG1;

analFolder = [testDir,'\AttDataInterp2\'];
frameNames = ls([analFolder,'\Frame*']);
%BAEdata = load([analFolder,frameNames(iImage,:)]);

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
%BAEdata = load([analFolder,frameNames(iImage,:)]);
%nLines = length(rayIdxs2);
%nKernels = size(IDF.Filt.F,2);

nEMLidx = 5;
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
                
            attInLiver = 2*liverThick * attLiver;
            attInLiverFilt = 2*liverThick * attLiverFilt;
            iDepth;
            
            corr = attLiverStruct.corr;
            corrFilt_COH = attLiverStruct.corr_filt;
        else
            attInLiver = nan; 
            attInLiverFilt = nan;
            corr = nan; 
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
 
            attLiverFilt = attLiverStruct.alpha_filt*1e6*vsxParams.Trans.frequency;
                
            attInLiverFilt_SNR = 2*liverThick * attLiverFilt; 
            iDepth;
            
            corr_SNR = attLiverStruct.corr;
            corrFilt_SNR = attLiverStruct.corr_filt;
        else
            attInLiver_SNR = nan; 
            attInLiverFilt_SNR = nan;
            corr_SNR = nan; 
            corrFilt_SNR = nan;
        end
    
        temp = [attInLiver,attInLiverFilt,attInLiverFilt_SNR;  
            corr,corrFilt_COH,corrFilt_SNR];
        
        if iEMLidx == 8 && iDepth == 2
            iEMLidx
        end

        for iCheck = 1:size(temp,2)
            if temp(1,iCheck)<0 || temp(1,iCheck)>10*(allDepths(iDepth)-liverStart*1e3).*1e-3*liverAtt
                temp(:,iCheck) = [nan,nan];
            end
        end

        
        %tempNew(iEMLidx,iDepth,:) = [attInLiver, attInLiverFilt,2* liverThick * liverAtt];

        attTest_DB = attSubCut + temp(1,1);
        attTest_DB_Filt= attSubCut + temp(1,2);
        
        attComp_Test =  10^(attTest_DB/10);
        attComp_Test_Filt =  10^(attTest_DB_Filt/10);
        
        allAtt_unseg(iDepth,:) = [attComp_Test,corr];
        allAtt_COH(iEMLidx,iDepth,:) = [attComp_Test_Filt,temp(2,2)];
    
        segBoolBIG_COH(iEMLidx,iDepth,:) = segBoolCluster_COH;

        attTest_DB_Filt_SNR = attSubCut + temp(1,3);
        
        attComp_Test_Filt =  10^(attTest_DB_Filt_SNR/10);
        
        allAtt_SNR(iEMLidx,iDepth,:) = [attComp_Test_Filt,temp(2,3)];
    
        segBoolBIG_SNR(iEMLidx,iDepth,:) = segBoolCluster_SNR;

    end
end

%% Attenuation plotting
close all

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


segPowCOH = nan([size(powf0BIG_COH)]);
segPowSNR = nan([size(powf0BIG_COH)]);

segPowCOH(segBoolBIG_COH) = powf0BIG_COH(segBoolBIG_COH); 
segPowSNR(segBoolBIG_SNR) = powf0BIG_SNR(segBoolBIG_SNR); 

bscCOH = zeros(nEMLidx,nDepth,2);
bscCOH(:,:,1) = mean(segPowCOH,3,'omitnan');
bscCOH(:,:,2) = std(segPowCOH,1,3,'omitnan');


bscSNR = zeros(nEMLidx,nDepth,2);
bscSNR(:,:,1) = mean(segPowSNR,3,'omitnan');
bscSNR(:,:,2) = std(segPowSNR,1,3,'omitnan');

powerUnseg = cell(1,length(allDepths));
powerSegCOHAll = cell(1,length(allDepths));
powerSegSNRAll = cell(1,length(allDepths));

nPossibleKernels = size(idxClustering(1:length(rayIdxs2),kWidth,oLap),2);    

for iDepth = 1:length(allDepths)

    segBoolDepthCOH = segBoolBIG_COH(:,iDepth,:);
    segBoolDepthSNR = segBoolBIG_SNR(:,iDepth,:);
    
    pctCOH = zeros(1,nEMLidx);
    pctSNR = zeros(1,nEMLidx);

    for iEML = 1:nEMLidx
        
        rayCoh = rayIdxs2(squeeze(segBoolDepthCOH(iEML,1,:)));
        kIdxsCOH = idxClustering(rayCoh,kWidth,oLap);
        nKernelsCOH = size(kIdxsCOH,2);

        raySNR = rayIdxs2(squeeze(segBoolDepthSNR(iEML,1,:)));
        kIdxsSNR = idxClustering(raySNR,kWidth,oLap);
        nKernelsSNR = size(kIdxsSNR,2);

        pctCOH(iEML) = 100*(1- nKernelsCOH/nPossibleKernels);
        pctSNR(iEML) = 100*(1- nKernelsSNR/nPossibleKernels);
    end

    bscUnseg = [mean(powf0BIG_unseg(:,iDepth)),std(powf0BIG_unseg(:,iDepth))];

    powerSegCOH = [[bscUnseg(1);bscCOH(:,iDepth,1)],[bscUnseg(2); bscCOH(:,iDepth,2)],nan(nEMLidx+1,1), [0;pctCOH'],nan(nEMLidx+1,1),nan(nEMLidx+1,1)]';
    powerSegSNR = [[bscUnseg(1);bscSNR(:,iDepth,1)],[bscUnseg(2); bscSNR(:,iDepth,2)],nan(nEMLidx+1,1), [0;pctSNR'],nan(nEMLidx+1,1),nan(nEMLidx+1,1)]';

    bscSegPlotterDual(powerSegCOH,powerSegSNR,1,0.1,mubsTH_STD,1)
    
    saveas(gcf,[imgDir,'/DoubleCorrBSC_',num2str(allDepths(iDepth)),'.fig'])
    exportgraphics(gcf,[imgDir,'/DoubleCorrBSC_',num2str(allDepths(iDepth)),'.pdf'])
    
    powerSegCOHAll{iDepth} = powerSegCOH;
    powerSegSNRAll{iDepth} = powerSegSNR;

end

save([imgDir,'/SegmResults.mat'],'powerSegSNRAll','powerSegCOHAll')
