
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

sumIdx = input('Sum Idx: ');


iImage = input('Which Image ? : ');


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
cohKlength_Length = 0.5*cohKlength*vsxParams.sPerWaveInit*vsxParams.lambda;


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

r1 = abs((skin.rho*skin.c - fat.rho*fat.c)/(skin.rho*skin.c + fat.rho*fat.c));
r2 = abs((muscle.rho*muscle.c - fat.rho*fat.c)/(muscle.rho*muscle.c + fat.rho*fat.c));

sAtt = skin.att(1)*8.6860000037 * vsxParams.Trans.frequency;
mAtt = muscle.att(1)*8.6860000037 * vsxParams.Trans.frequency;
fAtt = fat.att(1)*8.6860000037 * vsxParams.Trans.frequency;

attSubCut = skinThick*sAtt + fatThick*fAtt + muscleThick*mAtt;

attMeasures = load([testDir,'/AttDataTestFolderLiverDist/attFit',num2str(iImage),'.mat'])

attSpeckle  = [0.524 , 0.9].*vsxParams.Trans.frequency;
%attP        = [0.53  , 0.003].*vsxParams.Trans.frequency;


%ground truth BSC values: reference 
mu0_ref = 3.86e-4/(3^3.5); %(cm^{-1}sr^{-1}) 
bscSpeckleBf_ref = mu0_ref * 100 * vsxParams.Trans.frequency^3.5 ; %(m^{-1}sr^{-1})
bscSpeckleSTD_ref = [0.8 1.2].*mu0_ref

mubsTH = 0.1;
mubsTH_STD = [0.04 2.2];



%regions of image with full aperture in effect
xBool = 17:112;

rayIdxs = (1:128);    
rayIdxs2 = rayIdxs(xBool);  
EMLidx = 10;
kWidth = 5; 
oLap = 0.8;

allDepths = 15:5:50;

segBoolBIG1 = zeros(length(rayIdxs2),length(allDepths));
powf0_BIG = segBoolBIG1;
segBoolBIG2 = segBoolBIG1;

axIdxs_BIG = [];

iPlot=1;

figure;
swarmAx = gca;

allThValsCOH = zeros(2,length(allDepths));
allThValsSNR = zeros(2,length(allDepths));


for iDepths = 1:length(allDepths)

    speckleDir2 = [speckleDir,wName,'\Z',num2str(allDepths(iDepths)),'\'];  
    dataDir = [testDir,wName,'\Z',num2str(allDepths(iDepths)),'\'];
    

    dzT = allDepths(iDepths)*0.1 - edgeYVal*100; %cm 
    
    [~,attCoeffIDX] =  min(abs(dzT- attMeasures.xEst));
    
    liverThick = 2*(allDepths(iDepths)*1e-3 - (fatThick + muscleThick + skinThick))


    attLiver = -1.*(attMeasures.a0(attCoeffIDX)+attMeasures.alpha(attCoeffIDX)*vsxParams.Trans.frequency)
    
    attInLiver = liverThick * 100 * attLiver; 


    if attLiver <= 0 || isnan(attLiver)
        attInLiver = liverThick * 100 * 0.47 * vsxParams.Trans.frequency ;
    end
    
    attTest_DB = attSubCut + attInLiver;
    attComp_Test =  10^(attTest_DB/10);

    allAtt(iDepths) = attTest_DB;

    [~ , cohTestKernelIdx] = min(abs(allDepths(iDepths)*1e-3 - bfImgData.yVals));

    RMatTest = zeros(length(xBool),sumIdx);
    axIdxs = (cohTestKernelIdx- cohKlength/2):(cohTestKernelIdx+cohKlength/2-1);
    
    for iLine = 1:length(xBool)
       RMatTest(iLine,:) =  CoherenceAnalysisFN(squeeze(bfImgData.channelStack(xBool(iLine),axIdxs,1:sumIdx)));
    end
    
    axIdxs_BIG = [axIdxs_BIG,axIdxs];

    %% Coh analysis

    cohTest = sum(RMatTest(:,1:sumIdx),2);

    cInput = load([dataDir,'/Sum',num2str(sumIdx),'/COSIEinput',num2str(iImage),'.mat']);

    speckleCOSIE =  load([speckleDir2,'\COSIEoutput',adaptStr,num2str(cohKlength),'\COSIEoutput',num2str(sumIdx),'.mat']);
    cInput_Depth = load([dataDir,'/Sum',num2str(sumIdx),'/COSIEinput',num2str(iImage),'.mat']);

    powf0 = abs(cInput_Depth.spectAll(:,nF));
    specklePOWER_MEAN = mean(speckleCOSIE.powf0);
    
    convFactor = 1./specklePOWER_MEAN * bscSpeckleBf_ref *edgecorr*attComp_Test;

    bscEstimate = (powf0.*convFactor);
    powf0_BIG(:,iDepths) = bscEstimate;

   
    segBool1 = cohTest > speckleCOSIE.redEML(1,EMLidx) & cohTest < speckleCOSIE.redEML(2,EMLidx);    
    segBool1_cluster = ismember(rayIdxs2, unique(cell2mat(idxClustering(rayIdxs2(segBool1),kWidth,oLap))))';
    
    segBoolBIG1(:,iDepths) = segBool1_cluster;
    allThValsCOH(:,iDepths) = [speckleCOSIE.redEML(1,EMLidx),speckleCOSIE.redEML(2,EMLidx)];


   

    %% SNR analysis 

    qaSNRdata = load([testDir,wName,'\Z',num2str(allDepths(iDepths)),'\Sum',num2str(sumIdx),'\EnvStats',num2str(iImage),'.mat']);
    snrTest = qaSNRdata.envMean./qaSNRdata.envStd;
    
    %load COSIE data
    speckleSNRdata= load([speckleDir2 , 'envData_COSIE' ]) ;

    
    segBool2 = snrTest > speckleSNRdata.redEML(1,EMLidx) & snrTest < speckleSNRdata.redEML(2,EMLidx);
    segBool2_cluster = ismember(rayIdxs2, unique(cell2mat(idxClustering(rayIdxs2(segBool2),kWidth,oLap))))';
    segBoolBIG2(:,iDepths) = segBool2_cluster;
   
    
    allThValsSNR(:,iDepths) = [speckleSNRdata.redEML(1,EMLidx),speckleSNRdata.redEML(2,EMLidx)];
  
    %% 

    if 0

        hCohSpeckle = histogram(speckleCOSIE.thVector,'Normalization','probability');
        binCentreCoh = hCohSpeckle.BinEdges(2:end)-hCohSpeckle.BinWidth;
        pCoh = hCohSpeckle.Values;
        bscEstimate_weighted_COH = COVWeighting(cohTest,binCentreCoh,pCoh,bscEstimate,kWidth,oLap);
        
        weightingCOH(iDepths,:) = bscEstimate_weighted_COH(1:2,2);
        unseg(iDepths,:) = bscEstimate_weighted_COH(1:2,1);

        hEnvSpeckle = histogram(speckleSNRdata.snr,'Normalization','probability');
        binCentreEnv = hEnvSpeckle.BinEdges(2:end)-hEnvSpeckle.BinWidth;
        pEnv = hEnvSpeckle.Values;
        bscEstimate_weighted_ENV = COVWeighting(snrTest,binCentreEnv,pEnv,bscEstimate,kWidth,oLap);

        weightingENV(iDepths,:) = bscEstimate_weighted_ENV(1:2,2);
        
    end


    %% 

    

    if 1 %iDepths == 7 
        
        %bscSegPlotter(powerSeg,convFactor,bscSpeckleBf_test,bscSpeckleSTD_ref,dzT)
        dBSTD = 2.4;
        posErr = 0.1*(1+10^(dBSTD./10));
        negErr = -0.1*(1-10^(dBSTD./10));
       
        powerSegCOH= COVsegmentation_sK(cohTest,speckleCOSIE.EML,powf0,kWidth,oLap);
        %bscSegPlotter2(powerSegCOH,convFactor,0.1,[negErr,posErr],dzT)
       
        powerSegSNR = COVsegmentation_sK(snrTest,speckleSNRdata.EML,powf0,kWidth,oLap);
        %bscSegPlotter(powerSegSNR,convFactor,0.1,[negErr,posErr],dzT)

        bscSegPlotterDual(powerSegCOH,powerSegSNR,convFactor,mubsTH,mubsTH_STD,1)

        saveas(gcf,['Dual',num2str(dzT*10),'.fig'])
        saveas(gcf,['Dual',num2str(dzT*10),'.jpg'])
        exportgraphics(gcf,['Dual',num2str(dzT*10),'.pdf'], 'ContentType', 'vector');
          
        close all
    end

    %segDepthBoolIdxs = find(segBoolBIG1(:,iDepths));
    %segDepthBoolIdxs2 = find(~segBoolBIG1(:,iDepths));

    %fancySwarmPlotter(swarmAx,10*dzT,cohTest,segDepthBoolIdxs,speckleRej{iDepths},wireRej{iDepths},cystRej{iDepths})
    %fancySwarmPlotter(swarmAx,10*dzT,cohTest,segDepthBoolIdxs,segDepthBoolIdxs2,[],[])
    
    %fancySwarmPlotter(swarmAx,10*dzT,snrTest,segBool2_cluster,speckleRej{iDepths},wireRej{iDepths},cystRej{iDepths})

    if 0
        %f1 = figure;
        %a1 = axes;
        close all

        xPlot = xVals(xBool);
        xPlot = xPlot - mean(xPlot);
        xPlot = xPlot.*1e3;    
        %a1 = subplot(2,4,iDepths)
        a1 = axes
        plot(a1,xPlot,cohTest,'r-','LineWidth',2)
        hold on 
        plot(a1,[-16 16], [1 1].*speckleCOSIE.redEML(1,EMLidx),'k-.','LineWidth',2)
        plot(a1,[-16 16], [1 1].*speckleCOSIE.redEML(2,EMLidx),'k-.','LineWidth',2)
        xlim(a1,[-16 16])
        ylim(a1,[-10 30])
        xlabel(a1,'Lateral Position (mm)')
        ylabel(a1,'Coherence')
        set(a1,'FontSize',14)

        %saveas(gcf,['Coh',num2str(dzT*10),'mm.fig'])
        %exportgraphics(gcf, ['Coh',num2str(dzT*10),'mm.pdf'], 'ContentType', 'vector');
      
        close all
    end

     if 0
        %f1 = figure;
        close all

        a1 = axes;
        xPlot = xVals(xBool);
        xPlot = xPlot - mean(xPlot);
        xPlot = xPlot.*1e3;    
        %a1 = subplot(2,4,iDepths)
        plot(a1,xPlot,snrTest,'r-','LineWidth',2)
        hold on 
        plot(a1,[-16 16], [1 1].*speckleSNRdata.redEML(1,EMLidx),'k-.','LineWidth',2)
        plot(a1,[-16 16], [1 1].*speckleSNRdata.redEML(2,EMLidx),'k-.','LineWidth',2)
        xlim(a1,[-16 16])
        ylim(a1,[0 3.5])
        xlabel(a1,'Lateral Position (mm)')
        ylabel(a1,'Coherence')
        set(a1,'FontSize',14)

        saveas(gcf,['SNR',num2str(dzT*10),'mm.fig'])
        exportgraphics(gcf, ['SNR',num2str(dzT*10),'mm.pdf'], 'ContentType', 'vector');
      
        close all
    end



end

%plot(swarmAx,allDepths,allThValsCOH(1,:),'k-.','LineWidth',2)
%plot(swarmAx,allDepths,allThValsCOH(2,:),'k-.','LineWidth',2)

%plot(swarmAx,allDepths,allThValsSNR(1,:),'k-.','LineWidth',2)
%plot(swarmAx,allDepths,allThValsSNR(2,:),'k-.','LineWidth',2)

errorbar(15:5:50,unseg(:,1),unseg(:,2),'r.','MarkerSize',10)
hold on 
errorbar((15:5:50)-0.5,weightingCOH(:,1),weightingCOH(:,2),'k.','MarkerSize',10)
errorbar((15:5:50)+0.5 ,weightingENV(:,1),weightingENV(:,2),'b.','MarkerSize',10)
fillColor = [140, 222, 162]./256;
xShade = [5 5 60 60 ];
%yShade = bscSpeckleBf_test.*[0.8 1.2 1.2 0.8];
yShade = [mubsTH_STD(1) mubsTH_STD(2) mubsTH_STD(2) mubsTH_STD(1)];
%area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',bscSpeckleBf_test*.8,'ShowBaseLine','off');
area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',0.08,'ShowBaseLine','off');
plot([10 55],mubsTH.*[1 1],'-.','Color','b','LineWidth',2)
set(gca,'YScale','log')
set(gca,'FontSize',14)
xlabel('Depth (mm)')
xlim tight
ylabel('BSC (m^{-1}sr^{-1})')
saveas(gcf,'WeightingFigure.fig')
exportgraphics(gcf,['WeightingFigure.pdf'])



[ax1,ax2,cB] = bmCohCOSIEparImage_mDepth(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs_BIG,kLength_BSC_samples,powf0_BIG,segBoolBIG1,xBool,speckleCOSIE.pctSeg2(EMLidx),'Coherence');
xlim([-15 15])
ax1.XTick = [-15:5:15];
ax1.YTick = [2,5:10:55];
ax1.YLim = [2 55];
set(ax1,'FontSize',14)
set(ax1,'FontWeight','normal')
saveas(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage)])
saveas(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage),'.jpg'])
savefigPDF_Crop(gcf,[saveDir_SEG,'ParametricImage_COH',num2str(iImage)])


[ax3,ax4,cB2] = bmCohCOSIEparImage_mDepth(rayIdxs.*lWidth*1e3,yVals.*1e3,bfImgData,depthIdx,axIdxs_BIG,kLength_BSC_samples,powf0_BIG,segBoolBIG2,xBool,speckleSNRdata.pctSeg2(EMLidx),'Coherence');
xlim([-15 15])
ax3.XTick = [-15 :5:15];
ax3.YTick = [2 ,5:10:55];
ax3.YLim = [2 55];
set(ax3,'FontSize',18)
set(ax3,'FontWeight','normal')
saveas(gcf,[saveDir_SEG,'ParametricImage_SNR',num2str(iImage)])
saveas(gcf,[saveDir_SEG,'ParametricImage_SNR',num2str(iImage),'.jpg'])
savefigPDF_Crop(gcf,[saveDir_SEG,'ParametricImage_SNR',num2str(iImage)])





