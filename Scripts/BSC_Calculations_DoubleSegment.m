
clear 
close all 

path(path,'C:\Users\gwest\Documents\MATLAB\COSIE-Repo\Functions\')
path(path,'C:\Users\gwest\Documents\MATLAB\COSIE-Repo\Scripts\')
path(path,'C:\Users\gwest\Documents\MATLAB\AttenuationGUI\')


speckleDir = 'C:\Users\gwest\Documents\MATLAB\ElastPhtL74\Img1-4Dir\';

testDir = 'C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\';

%load verasonics params2
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);



allDepths = 20:5:55;

nImages = 7;
nDepths = length(allDepths);
nEML = 10;

powerUnseg = zeros(nImages,nDepths,2);
powerSegCOH = zeros(nImages,nDepths,nEML,2);
powerSegSNR = powerSegCOH;

for iImage = 1:nImages

   
    imgDir = ['C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\SegmResults\'];
    s1 = load([imgDir,'\SegmResults',num2str(iImage),'.mat']);
    
    for iDepth = 1:nDepths

        temp1 = s1.powerSegCOHAll{iDepth};
        temp2 = s1.powerSegSNRAll{iDepth};

        powerUnseg(iImage,iDepth,:) = [temp1(1,1) temp1(2,1)];

        powerSegCOH(iImage,iDepth,:,1) = temp1(1,2:(nEML+1));
        powerSegCOH(iImage,iDepth,:,2) = temp1(2,2:(nEML+1));

        powerSegSNR(iImage,iDepth,:,1) = temp2(1,2:(nEML+1));
        powerSegSNR(iImage,iDepth,:,2) = temp2(2,2:(nEML+1));
     
    end
end


%% all image mean

powerUnsegMean = squeeze(mean(powerUnseg(:,:,1),1,'omitnan'));
powerUnsegStd  = squeeze(1/nImages.*sqrt(sum(powerUnseg(:,:,2).^2,1,'omitnan')));

powerSegMean = squeeze(mean(powerSegCOH(:,:,:,1),1,'omitnan'));
powerSegStd = squeeze(1/nImages.*sqrt(sum(powerSegCOH(:,:,:,2).^2,1,'omitnan')));


load('C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\SegmResults\CMap.mat')



fillColor = [140, 222, 162]./256;

mubsTH = 0.1;
mubsTH_STD = [0.04 2.2];
neg = mubsTH_STD(1);
pos = mubsTH_STD(2);




% Coh segmentation, depth plot
offSet = 1.5;
[powerUnsegMean1,powerUnsegStd1,powerSegMean_COH,powerSegStd_COH] = doubleSegmentationAveragePlot(powerUnseg,powerSegCOH,linspace(-offSet,offSet,nEML+1),allDepths,mubsTH, mubsTH_STD,1,1:8,1:7,cmap);
title('Coherence Segm')
saveas(gcf,[imgDir,'/COHdepthSeg.jpg' ])


% Coherence Segmentation, image plot
offSet = 0.3; 
doubleSegmentationAveragePlot(powerUnseg,powerSegCOH,linspace(-offSet,offSet,nEML+1),1:7,mubsTH, mubsTH_STD,2,1:7,1:7,cmap);
colormap([0,0,0;cmap])
cB = colorbar;
set(cB,'XTick',(1:2:(2*(nEML+1)))./(2*(nEML+1)))
set(cB,'TickLabels',0:nEML)
title('Coherence Segm')
saveas(gcf,[imgDir,'/COHimageSeg.jpg' ])

% SNR segmentation, depth plot 
doubleSegmentationAveragePlot(powerUnseg,powerSegSNR,linspace(-offSet,offSet,nEML+1),allDepths,mubsTH, mubsTH_STD,1,1:8,1:7,cmap);
title('SNR Segm')

saveas(gcf,[imgDir,'/SNRdepthSeg.jpg' ])


% SNR segmentation, depth plot
doubleSegmentationAveragePlot(powerUnseg,powerSegSNR,linspace(-offSet,offSet,nEML+1),1:7,mubsTH, mubsTH_STD,2,1:7,1:7,cmap);
colormap([0,0,0;cmap])
cB = colorbar;
set(cB,'XTick',(1:2:(2*(nEML+1)))./(2*(nEML+1)))
set(cB,'TickLabels',0:nEML)
title('SNR Segm')
saveas(gcf,[imgDir,'/SNRimageSeg.jpg' ])



%% Full averaging 


nEML = 10;

powerUnseg = zeros(nImages,nDepths,2);
powerSegCOH = zeros(nImages,nDepths,nEML,2);
powerSegSNR = powerSegCOH;

for iImage = 1:nImages

   
    %imgDir = ['C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\SegmResults\'];
    s1 = load([imgDir,'\SegmResults',num2str(iImage),'.mat']);
    
    for iDepth = 1:nDepths

        temp1 = s1.powerSegCOHAll{iDepth};
        temp2 = s1.powerSegSNRAll{iDepth};

        powerUnseg(iImage,iDepth,:) = [temp1(1,1) temp1(2,1)];

        powerSegCOH(iImage,iDepth,:,1) = temp1(1,2:(nEML+1));
        powerSegCOH(iImage,iDepth,:,2) = temp1(2,2:(nEML+1));

        powerSegSNR(iImage,iDepth,:,1) = temp2(1,2:(nEML+1));
        powerSegSNR(iImage,iDepth,:,2) = temp2(2,2:(nEML+1));
     
    end
end

powerUnsegMean = squeeze(mean(powerUnseg(:,:,1),1,'omitnan'));
powerUnsegStd  = squeeze(1/nImages.*sqrt(sum(powerUnseg(:,:,2).^2,1,'omitnan')));

powerSegMean_COH = squeeze(mean(powerSegCOH(:,:,:,1),1,'omitnan'));
powerSegStd_COH = squeeze(1/nImages.*sqrt(sum(powerSegCOH(:,:,:,2).^2,1,'omitnan')));

powerSegMean_SNR = squeeze(mean(powerSegSNR(:,:,:,1),1,'omitnan'));
powerSegStd_SNR = squeeze(1/nImages.*sqrt(sum(powerSegSNR(:,:,:,2).^2,1,'omitnan')));



nNormUnseg = sum(~isnan(powerUnsegMean1(:,:,1)));
powerUnsegMean2 = mean(powerUnsegMean1(:,:,1),'omitnan');
powerUnsegStd2  = squeeze(1./nNormUnseg.*sqrt(sum(powerUnsegStd1.^2,'omitnan')));

   
nNorm = sum(~isnan(powerSegMean_COH(:,:,1)),1);
powerSegMean_COH2 = squeeze(mean(powerSegMean_COH,1,'omitnan'));
powerSegStd_COH2 = squeeze(1./nNorm.*sqrt(sum(powerSegStd_COH.^2,1,'omitnan')));

   
nNorm = sum(~isnan(powerSegMean_SNR(:,:,1)),1);
powerSegMean_SNR2 = squeeze(mean(powerSegMean_SNR,1,'omitnan'));
powerSegStd_SNR2 = squeeze(1./nNorm.*sqrt(sum(powerSegStd_SNR.^2,1,'omitnan')));


figure
errorbar(0,powerUnsegMean2,powerUnsegStd2,'ro','MarkerFaceColor','r')
hold on 
errorbar(1:10,powerSegMean_COH2,powerSegStd_COH2,'k-o','MarkerFaceColor','k')
errorbar(1:10,powerSegMean_SNR2,powerSegStd_SNR2,'b-o','MarkerFaceColor','b')
xlim padded

xlabel('EML index ')
ylabel('BSC (m^{-1}rad^{-1})')
set(gca,'YScale','log')
grid on 

saveas(gcf,[imgDir,'/EMLSeg.jpg' ])

