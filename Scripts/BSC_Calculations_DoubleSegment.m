
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

allDepths = 15:5:50;

nImages = 6;
nDepths = length(allDepths);
nEML = 5;

powerUnseg = zeros(nImages,nDepths,2);
powerSegCOH = zeros(nImages,nDepths,nEML,2);
powerSegSNR = powerSegCOH;

for iImage = 1:nImages

    imgDir = ['C:\Users\gwest\Documents\COSIE paper 1\EmmaLiver\Img',num2str(iImage),'\AttSeg\'];
    s1 = load([imgDir,'\SegmResults.mat']);
    
    for iDepth = 1:8

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

cmap = [
    0.1216 0.4667 0.7059;  % Strong Blue
    0.8510 0.3725 0.0078;  % Vivid Orange
    0.5569 0.2667 0.6784; 
    0.7922 0.1176 0.5686;  % Magenta
    0.4980 0.4980 0.4980;  % Dark Gray
];

fillColor = [140, 222, 162]./256;

mubsTH = 0.1;
mubsTH_STD = [0.04 2.2];
neg = mubsTH_STD(1);
pos = mubsTH_STD(2);

offSet = 1.5;
doubleSegmentationAveragePlot(powerUnseg,powerSegCOH,linspace(-offSet,offSet,nEML+1),allDepths,mubsTH, mubsTH_STD,1,1:8,1:6,cmap)
doubleSegmentationAveragePlot(powerUnseg,powerSegSNR,linspace(-offSet,offSet,nEML+1),allDepths,mubsTH, mubsTH_STD,1,1:8,1:6,cmap)

offSet = 0.3; 
doubleSegmentationAveragePlot(powerUnseg,powerSegCOH,linspace(-offSet,offSet,nEML+1),1:6,mubsTH, mubsTH_STD,2,1:6,1:6,cmap)
colormap([0,0,0;cmap])
cB = colorbar;
set(cB,'XTick',(1:2:(2*(nEML+1)))./(2*(nEML+1)))
set(cB,'TickLabels',0:nEML)


doubleSegmentationAveragePlot(powerUnseg,powerSegSNR,linspace(-offSet,offSet,nEML+1),1:6,mubsTH, mubsTH_STD,2,1:6,1:6,cmap)
colormap([0,0,0;cmap])
cB = colorbar;
set(cB,'XTick',(1:2:(2*(nEML+1)))./(2*(nEML+1)))
set(cB,'TickLabels',0:nEML)



    

