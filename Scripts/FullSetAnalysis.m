clear 
close all

topDir = 'C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\';
resultsDir = 'C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\SegmResults\';

cmapStruct= load([resultsDir,'\CMap']);

vsxParams = load([topDir,'VSXoutput.mat']);

fNames = ls([resultsDir,'SegmResults*']);

mubsTH = 0.1;
mubsTH_STD = [0.04 2.2];

nImages = size(fNames,1);

dataAll = cell(1,nImages);

for iImage = 1 :nImages 

    dataAll{iImage} = load([resultsDir,fNames(iImage,:)]);

end 

%% Average attenuation plots 

nEML = size(dataAll{1}.allAtt_SNR,1);


allAttCoh = zeros(nImages,size(dataAll{1}.allAtt_COH,1),size(dataAll{1}.allAtt_COH,2),2);
allAttSNR = allAttCoh;
allAttUnseg = zeros(nImages,size(dataAll{1}.allAtt_COH,2) ,2);

for iImage = 1:nImages
    %for iEML = 1:10
        
        allAttCoh(iImage,:,:,:) = dataAll{iImage}.allAtt_COH;
        allAttSNR(iImage,:,:,:) = dataAll{iImage}.allAtt_SNR;
        allAttUnseg(iImage,:,:) = dataAll{iImage}.allAtt_unseg;
    
    %end
end


%% Mean attenuation by image dimension 


meanAttCoh = squeeze(mean(allAttCoh,1,'omitnan'));
meanAttSNR = squeeze(mean(allAttSNR,1,'omitnan'));
meanAttUnseg = squeeze(mean(allAttUnseg,1,'omitnan'));


stdAttCoh = squeeze(std(allAttCoh,1,1,'omitnan'));
stdAttSNR = squeeze(std(allAttSNR,1,1,'omitnan'));
stdAttUnseg = squeeze(std(allAttUnseg,1,1,'omitnan'));


[f1,f2] = AverageAttenuationPlotter(dataAll,meanAttCoh,meanAttUnseg,cmapStruct.cmap);

saveas(f1,[resultsDir,'\MeanAttValuesCOH.fig' ])
saveas(f2,[resultsDir,'\ResidValuesCOH.fig' ])


[f3,f4] = AverageAttenuationPlotter(dataAll,meanAttSNR,meanAttUnseg,cmapStruct.cmap);

saveas(f3,[resultsDir,'\MeanAttValuesSNR.fig' ])
saveas(f4,[resultsDir,'\ResidValuesSNR.fig' ])

close all


%% Mean BSC by image dimension



for iDepth = 1 : length(dataAll{1}.allDepths)   

    bscSegCOH = zeros(size(dataAll{1}.powerSegCOHAll{1},1),size(dataAll{1}.powerSegCOHAll{1},2),nImages);    
    bscSegSNR = bscSegCOH;

    for iImage = 1:nImages
        
        bscSegCOH(:,:,iImage) = dataAll{iImage}.powerSegCOHAll{iDepth};    
        bscSegSNR(:,:,iImage) = dataAll{iImage}.powerSegSNRAll{iDepth};    
            
    end

    bscSegCOHmean = mean(bscSegCOH,3,'omitnan');
    bscSegSNRmean = mean(bscSegSNR,3,'omitnan');
    
    bscSegPlotterDual(bscSegCOHmean,bscSegSNRmean,1,mubsTH,mubsTH_STD,1)
    title(['Z = ',num2str(dataAll{1}.allDepths(iDepth))])
    saveas(gcf,[resultsDir,'\MeanBSC',num2str(iDepth),'.fig'])

    unseg(iDepth,:) = bscSegCOHmean(1:2,1);
 
end