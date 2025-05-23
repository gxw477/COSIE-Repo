
clear 
close all 


path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE\COSIE-Repo\Functions\')
path(path,'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE\COSIE-Repo\Scripts\')


speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\EmmaLiver\EmmaLiver_NHV_NTGC\';
planeDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Pref\';

%load verasonics param's2
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

sumIdx = 33;

%Calculate kernel idxs for analysis
kLength_BSC_samples = 120;
cohKlength = kLength_BSC_samples;
cohKlength_Length = 0.5*cohKlength*vsxParams.sPerWaveInit*vsxParams.lambda;
kWidth = 5; 
oLap = 0.8;

samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample  + 1;
yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;

lambda = 1540/vsxParams.Trans.frequency*1e-6;
lWidth = lambda*(vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1)) ;
xVals = (1:128).*lWidth;


%spectral calculations
fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;
fVals = (0:kLength_BSC_samples-1).*fs/(kLength_BSC_samples);
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
end

%% Ground truth values 

%Emma Liver 

attData = load([testDir,'AttData\attFit',num2str(1),'.mat']);
attTest = vsxParams.Trans.frequency.*attData.alpha+attData.a0;



%Reference phantom
attSpeckle  = [0.524 , 0.9].*vsxParams.Trans.frequency;

%ground truth BSC values: reference 
mu0_ref = 3.86e-4/(3^3.5); %(cm^{-1}sr^{-1}) 
bscSpeckleBf_ref = mu0_ref * 100 * vsxParams.Trans.frequency^3.5 ; %(m^{-1}sr^{-1})
bscSpeckleSTD_ref = 0.2*bscSpeckleBf_ref;

%ground truth BSC values: test
bscSpeckleBf_test = 0.1; %(m^{-1}sr^{-1})
bscSpeckleSTD_test = 0.2*bscSpeckleBf_test;

%% 

%regions of image with full aperture in effect
xBool = 17:112;
rayIdxs = (1:128);    
rayIdxs2 = rayIdxs(xBool);  
EMLidx = 10;
allDepths = 15:5:50;
nDepths = length(allDepths);

fNames = ls([testDir,'BFimg*']);

nImages = size(fNames,1);



bscResults = cell(nImages,nDepths);


for iImage = 1:nImages

    iImage

    %load test data
    bfImgData = load([testDir,fNames(iImage,:)]);
    
    adaptBool = 1;%input('Adaptive grid sizing? : ');
    
    if adaptBool 
        adaptStr = '_adaptive';
    else
        adaptStr = '';
    end
    
    wOption = 2;%input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch \n 4 for Hanning: \n ');
    
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
    
    %test data interface reflectivity
    bm = bmode(bfImgData.iq',80);  
    [~ , edge] = (max(bm,[],1));
    edgeYVal = 0;%bfImgData.yVals(round(mean(edge)));
   

    for iDepths = 1:nDepths

        speckleDir2 = [speckleDir,wName,'\Z',num2str(allDepths(iDepths)),'\'];  
        dataDir = [testDir,wName,'\Z',num2str(allDepths(iDepths)),'\'];
        
    
        dzT = allDepths(iDepths)*0.1 - edgeYVal*100; %cm 
        
        [~,idx] = min(abs(attData.xEst-dzT));

        attTest_DB = 2*(dzT)*attTest(idx);
        attComp_Test =   10^(attTest_DB/10);
    
        [~ , cohTestKernelIdx] = min(abs(allDepths(iDepths)*1e-3 - bfImgData.yVals));
    
        RMatTest = zeros(length(xBool),sumIdx);
        axIdxs = (cohTestKernelIdx- cohKlength/2):(cohTestKernelIdx+cohKlength/2-1);
        
        for iLine = 1:length(xBool)
           RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(xBool(iLine),axIdxs,:)));
        end 
        
        cInput_Depth = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);
        powf0 = abs(cInput_Depth.spectAll(:,nF));
        
        speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(cohKlength),'\COSIEoutput',num2str(sumIdx),'.mat']);
        specklePOWER_MEAN = mean(speckleCOSIE.powf0);
        convFactor = .1./specklePOWER_MEAN * bscSpeckleBf_ref *edgecorr*attComp_Test;
    
        bscEstimate = (powf0.*convFactor);
        powf0_BIG(:,iDepths) = bscEstimate;
        
        
        %% Coh analysis
        
        cohTest = sum(RMatTest(:,1:sumIdx),2);    
        cInput = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);
        cohSeg = COVsegmentation_sK(cohTest,speckleCOSIE.EML,bscEstimate,kWidth,oLap);


        %% SNR analysis 
        
        %load COSIE data
        speckleSNRdata= load([speckleDir2 , 'envData_COSIE' ]) ;
        qaSNRdata = load([testDir,wName,'\Z',num2str(allDepths(iDepths)),'\EnvStats',num2str(iImage),'.mat']);
        snrTest = qaSNRdata.envMean./qaSNRdata.envStd;
        
        snrSeg = COVsegmentation_sK(snrTest,speckleSNRdata.EML,bscEstimate,kWidth,oLap);

        
      
        %% Weighting analysis 

        hCohSpeckle = histogram(speckleCOSIE.thVector,'Normalization','probability');
        binCentreCoh = hCohSpeckle.BinEdges(2:end)-hCohSpeckle.BinWidth;
        pCoh = hCohSpeckle.Values;
        bscEstimate_weighted_COH = COVWeighting(cohTest,binCentreCoh,pCoh,bscEstimate,kWidth,oLap);

        close all
  
        hEnvSpeckle = histogram(speckleSNRdata.snr,'Normalization','probability');
        binCentreEnv = hEnvSpeckle.BinEdges(2:end)-hEnvSpeckle.BinWidth;
        pEnv = hEnvSpeckle.Values;
        bscEstimate_weighted_ENV = COVWeighting(snrTest,binCentreEnv,pEnv,bscEstimate,kWidth,oLap);

        close all

        %% Store results
        
        bscResultsIdepth.coherenceCOSIE = cohSeg;
        bscResultsIdepth.snrCOSIE = snrSeg;
        bscResultsIdepth.coherenceWEIGHT = bscEstimate_weighted_COH;
        bscResultsIdepth.snrWEIGHT = bscEstimate_weighted_ENV;
        

        bscResults{iImage,iDepths} = bscResultsIdepth;


    end

end



%% depth 

%already set but leaving here to allow change
EMLidx = 10; 

cohCOSIEdepth = zeros(nDepths,2);
snrCOSIEdepth = zeros(nDepths,2);
cohWEIGHTdepth = zeros(nDepths,2);
snrWEIGHTdepth =  zeros(nDepths,2);
unseg = zeros(nDepths,2);

for iDepth = 1:nDepths
    
    cohCOSIEtemp = zeros(nImages,2);
    snrCOSIEtemp = cohCOSIEtemp;
    cohWEIGHTtemp = snrCOSIEtemp;
    snrWEIGHTtemp = cohWEIGHTtemp;
    unsegTemp = snrWEIGHTtemp;
    
    for iImage = 1:nImages

        cohCOSIEtemp(iImage,:)  = (bscResults{iImage,iDepth}.coherenceCOSIE([1,2],EMLidx));
        snrCOSIEtemp(iImage,:)  = (bscResults{iImage,iDepth}.snrCOSIE([1,2],EMLidx));
        cohWEIGHTtemp(iImage,:) = (bscResults{iImage,iDepth}.coherenceWEIGHT([1,2],2));
        snrWEIGHTtemp(iImage,:) = (bscResults{iImage,iDepth}.snrWEIGHT([1,2],2));
        unsegTemp(iImage,:)     = (bscResults{iImage,iDepth}.snrWEIGHT([1,2],1));
        
    end
   
    cohCOSIEdepth(iDepth,:) = [mean(cohCOSIEtemp(:,1),'omitnan'),mean(cohCOSIEtemp(:,2),'omitnan')];
    snrCOSIEdepth(iDepth,:) = [mean(snrCOSIEtemp(:,1),'omitnan'),mean(snrCOSIEtemp(:,2),'omitnan')];
    cohWEIGHTdepth(iDepth,:) = [mean(cohWEIGHTtemp(:,1),'omitnan'),mean(cohWEIGHTtemp(:,2),'omitnan')];
    snrWEIGHTdepth(iDepth,:) = [mean(snrWEIGHTtemp(:,1),'omitnan'),mean(snrWEIGHTtemp(:,2),'omitnan')];
    unseg(iDepth,:) = [mean(unsegTemp(:,1),'omitnan'),mean(unsegTemp(:,2),'omitnan')];
    
end


figure 
errorbar((15:5:50)-.5,cohCOSIEdepth(:,1),cohCOSIEdepth(:,2),'k.','MarkerSize',10)
hold on 
errorbar((15:5:50)+.2,snrCOSIEdepth(:,1),snrCOSIEdepth(:,2),'r.','MarkerSize',10)
%errorbar((15:5:50)+.2,cohWEIGHTdepth(:,1),cohWEIGHTdepth(:,2),'k*','MarkerSize',10)
%errorbar((15:5:50)+.5,snrWEIGHTdepth(:,1),snrWEIGHTdepth(:,2),'r*','MarkerSize',10)

fillColor = [16,158,99]./256;
errorbar((15:5:50),unseg(:,1),unseg(:,2),'o','MarkerSize',2,'Color',fillColor,'MarkerFaceColor',fillColor)
set(gca,'YScale','log')


figure 
plot((15:5:50)-.5,cohCOSIEdepth(:,2),'k.','MarkerSize',10)
hold on 
plot((15:5:50)+.2,snrCOSIEdepth(:,2),'r.','MarkerSize',10)
plot((15:5:50)+.2,cohWEIGHTdepth(:,2),'k*','MarkerSize',10)
plot((15:5:50)+.5,snrWEIGHTdepth(:,2),'r*','MarkerSize',10)

fillColor = [16,158,99]./256;
plot((15:5:50),unseg(:,2),'o','MarkerSize',2,'Color',fillColor,'MarkerFaceColor',fillColor)
set(gca,'YScale','log')





