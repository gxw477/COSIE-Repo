
clear 
close all 


path(path,'C:\Users\gwest\Documents\MATLAB\COSIE-Repo\Functions\')
path(path,'C:\Users\gwest\Documents\MATLAB\COSIE-Repo\Scripts\')
path(path,'C:\Users\gwest\Documents\MATLAB\AttenuationGUI\')



speckleDir = 'C:\Users\gwest\Documents\MATLAB\ElastPhtL74\Img1-4Dir\QUAD\';
testDir = 'C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\';
%planeDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Pref\';

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

r1 = abs((skin.rho*skin.c - fat.rho*fat.c)/(skin.rho*skin.c + fat.rho*fat.c));
r2 = abs((muscle.rho*muscle.c - fat.rho*fat.c)/(muscle.rho*muscle.c + fat.rho*fat.c));

sAtt = skin.att(1)*8.6860000037 * vsxParams.Trans.frequency;
mAtt = muscle.att(1)*8.6860000037 * vsxParams.Trans.frequency;
fAtt = fat.att(1)*8.6860000037 * vsxParams.Trans.frequency;

attSubCut = skinThick*sAtt + fatThick*fAtt + muscleThick*mAtt;


attSpeckle  = [0.524 , 0.9].*vsxParams.Trans.frequency;
%attP        = [0.53  , 0.003].*vsxParams.Trans.frequency;


%ground truth BSC values: reference 
mu0_ref = 3.86e-4/(3^3.5); %(cm^{-1}sr^{-1}) 
bscSpeckleBf_ref = mu0_ref * 100 * vsxParams.Trans.frequency^3.5 ; %(m^{-1}sr^{-1})
bscSpeckleSTD_ref = 0.2*bscSpeckleBf_ref;


%% 

%regions of image with full aperture in effect
xBool = 17:112;
rayIdxs = (1:128);    
rayIdxs2 = rayIdxs(xBool);  
EMLidx = 5;
allDepths = 20:5:55;
nDepths = length(allDepths);

fNames = ls([testDir,'BFimg*']);

nImages = size(fNames,1);

bscResults = cell(nImages,nDepths);

axIdxs_BIG = [];


for iImage = 1:nImages

    iImage

    %load test data
    bfImgData = load([testDir,fNames(iImage,:)]);
    
    %attMeasures = load([testDir,'/AttData/Frame',num2str(iImage),'.mat'])
    attMeasures = load([testDir,'\AttDataTestFolderLiverDist\attFit',num2str(iImage),'.mat']);

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
        dataDir = [testDir,wName,'\Z',num2str(allDepths(iDepths)),'\COSIEoutput_adaptive120\Sum33\'];
        
    
        dzT = allDepths(iDepths)*0.1 - edgeYVal*100; %cm 
        
        [~,attCoeffIDX] =  min(abs(dzT- attMeasures.xEst));
        
        liverThick = allDepths(iDepths)*1e-3 - (fatThick + muscleThick + skinThick)
        
        
        attLiver = -1.*(attMeasures.a0(attCoeffIDX)+attMeasures.alpha(attCoeffIDX)*vsxParams.Trans.frequency)
        
        attInLiver = liverThick * 100 * attLiver; 
        
        
        if attLiver < 0 || isnan(attLiver)
            attInLiver = liverThick * 100 * 0.47 * vsxParams.Trans.frequency ;
        end
        
        attTest_DB = attSubCut + attInLiver;
        attComp_Test =  10^(attTest_DB/10);
        
        allAtt(iImage,iDepths) = attTest_DB;
        
        [~ , cohTestKernelIdx] = min(abs(allDepths(iDepths)*1e-3 - bfImgData.yVals));
        
        RMatTest = zeros(length(xBool),sumIdx);
        axIdxs = (cohTestKernelIdx- cohKlength/2):(cohTestKernelIdx+cohKlength/2-1);
        
        for iLine = 1:length(xBool)
           RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(xBool(iLine),axIdxs,:)));
        end 
        axIdxs_BIG = [axIdxs_BIG,axIdxs];

        %% Coh Analysis
        
        cohTest = sum(RMatTest(:,1:sumIdx),2);

        cInput_Depth = load([dataDir,'COSIEinput',num2str(iImage),'.mat']);
        powf0 = abs(cInput_Depth.spectAll(:,nF));
        
        speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(cohKlength),'\COSIEoutput',num2str(sumIdx),'.mat']);
        specklePOWER_MEAN = mean(speckleCOSIE.powf0);
        convFactor = 1./specklePOWER_MEAN * bscSpeckleBf_ref *edgecorr*attComp_Test;
    
        bscEstimate = (powf0.*convFactor);
        powf0_BIG(:,iDepths) = bscEstimate;
        
        

        %% SNR analysis 
        
        qaSNRdata = load([testDir,wName,'\Z',num2str(allDepths(iDepths)),'\COSIEoutput_adaptive',num2str(cohKlength),'\Sum',num2str(sumIdx),'\EnvStats',num2str(iImage),'.mat']);
        snrTest = qaSNRdata.envMean./qaSNRdata.envStd;
        
        %load COSIE data
        speckleSNRdata= load([speckleDir2 , 'envData_COSIE' ]) ;
    
        
        segBool2 = snrTest > speckleSNRdata.redEML(1,EMLidx) & snrTest < speckleSNRdata.redEML(2,EMLidx);
        segBool2_cluster = ismember(rayIdxs2, unique(cell2mat(idxClustering(rayIdxs2(segBool2),kWidth,oLap))))';
        segBoolBIG2(:,iDepths) = segBool2_cluster;
       
        allThValsSNR(:,iDepths) = [speckleSNRdata.redEML(1,EMLidx),speckleSNRdata.redEML(2,EMLidx)];
        
      
        %% Weighting analysis 

        hCohSpeckle = histogram(speckleCOSIE.thVector,'Normalization','probability');
        binCentreCoh = hCohSpeckle.BinEdges(2:end)-hCohSpeckle.BinWidth;
        pCoh = hCohSpeckle.Values;
        bscEstimate_weighted_COH = COVWeighting(cohTest,binCentreCoh,pCoh,bscEstimate,kWidth,oLap);

        close(gcf)
  
        hEnvSpeckle = histogram(speckleSNRdata.snr,'Normalization','probability');
        binCentreEnv = hEnvSpeckle.BinEdges(2:end)-hEnvSpeckle.BinWidth;
        pEnv = hEnvSpeckle.Values;
        bscEstimate_weighted_ENV = COVWeighting(snrTest,binCentreEnv,pEnv,bscEstimate,kWidth,oLap);

        close(gcf)

        %% Store results

        cohSeg = COVsegmentation_sK(cohTest,speckleCOSIE.redEML,powf0,kWidth,oLap);
        snrSeg = COVsegmentation_sK(snrTest,speckleSNRdata.redEML,powf0,kWidth,oLap);
      
        bscResultsIdepth.coherenceCOSIE = cohSeg;
        bscResultsIdepth.snrCOSIE = snrSeg;
        bscResultsIdepth.coherenceWEIGHT = bscEstimate_weighted_COH;
        bscResultsIdepth.snrWEIGHT = bscEstimate_weighted_ENV;
        

        bscResults{iImage,iDepths} = bscResultsIdepth;


    end

    %plot(allAtt)
    %hold on 
end

        

%% depth 

%already set but leaving here to allow change
EMLidx = 5 

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
        unsegTemp(iImage,:)     = (bscResults{iImage,iDepth}.coherenceCOSIE([1,2],1));
        
    end
   
    cohCOSIEdepth(iDepth,:) = [mean(cohCOSIEtemp(:,1),'omitnan'),mean(cohCOSIEtemp(:,2),'omitnan')].*convFactor;
    snrCOSIEdepth(iDepth,:) = [mean(snrCOSIEtemp(:,1),'omitnan'),mean(snrCOSIEtemp(:,2),'omitnan')].*convFactor;
    cohWEIGHTdepth(iDepth,:) = [mean(cohWEIGHTtemp(:,1),'omitnan'),mean(cohWEIGHTtemp(:,2),'omitnan')];
    snrWEIGHTdepth(iDepth,:) = [mean(snrWEIGHTtemp(:,1),'omitnan'),mean(snrWEIGHTtemp(:,2),'omitnan')];
    unseg(iDepth,:) = [mean(unsegTemp(:,1),'omitnan'),mean(unsegTemp(:,2),'omitnan')].*convFactor;
    
end


mubsTH = 0.1;
dBSTD = 2.4;
posErr = mubsTH*(1+10^(dBSTD./10));
negErr = -mubsTH*(1-10^(dBSTD./10));

bscSpeckleSTD_ref = [0.04 2.2]

fillColor = [140, 222, 162]./256;
xShade = [-15 -15 110 110 ];
yShade = mubsTH + [- negErr posErr posErr -negErr ];
    


%%

figure 
area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',mubsTH-negErr,'ShowBaseLine','off');    
hold on 
plot([-15 110],mubsTH.*[1,1],'-.','Color',[255 171 0]./256,'LineWidth',2)
errorbar(allDepths,unseg(:,1),unseg(:,2),'r.')
errorbar(allDepths-1,cohCOSIEdepth(:,1),cohCOSIEdepth(:,2),'k.')
errorbar(allDepths-0.5,snrCOSIEdepth(:,1),snrCOSIEdepth(:,2),'b.')
errorbar(allDepths+0.5,cohWEIGHTdepth(:,1),cohWEIGHTdepth(:,2),'kx','MarkerFaceColor','k','MarkerSize',5)
errorbar(allDepths+1,snrWEIGHTdepth(:,1),snrWEIGHTdepth(:,2),'bx','MarkerFaceColor','b','MarkerSize',5)
set(gca,'YScale','log')
xlim([10 55])

xlabel('Depth (mm)')
ylabel('log_{10}(BSC)')
set(gca,'FontSize',20)
ylim([1e-5-1e-6 1e1])
%yticklabels({'-5','-4','-3','-2','-1','0','1'})
yticks([10.^[-5: 1]])

%% Comparison between SNR COSIE and coherent COSIE

compMeasure = [cohCOSIEdepth(:,1),snrCOSIEdepth(:,1)];

plot(cohCOSIEdepth(:,1),snrCOSIEdepth(:,1),'ko','MarkerFaceColor','k')
hold on 

xFit= cohCOSIEdepth(2:end,1);
yFit = snrCOSIEdepth(2:end,1);

[a,b,corr,err,res] = lsqn(xFit(2:end),yFit(2:end),[0 0 0]);

xPlot = linspace(0,max(compMeasure(:)),100);
yPlot = xPlot.*a + b;

plot(xPlot,yPlot,'k-.')


%% 

unsegCOV = unseg(:,2)./unseg(:,1);
cohCOV = cohCOSIEdepth(:,2)./cohCOSIEdepth(:,1);
snrCOV = snrCOSIEdepth(:,2)./snrCOSIEdepth(:,1);


figure 
plot(allDepths,unsegCOV,'ro','MarkerFaceColor','r')
hold on 
plot(allDepths-1,cohCOV,'ko','MarkerFaceColor','k')
plot(allDepths-0.5,snrCOV,'bo','MarkerFaceColor','b')
xlabel('Depth (mm)')
ylabel('C.O.V.')
set(gca,'FontSize',20)
xlim([10 55])
