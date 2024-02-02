


close all 
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

zSelect = input('Depth of interest (mm) : ');

speckleDir = [speckleDir,'Z',num2str(zSelect),'\'];
topDir2 = [topDir,'Z',num2str(zSelect),'\'];
testData = load([topDir2,'COSIEinput',num2str(iImage),'.mat']);


tgcBool = isequal(vsxParams.TGC,vsxParams2.TGC);

if ~tgcBool
    error('TGC altered between Speckle and test image')
end


lWidth = (vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1))*vsxParams.lambda;
xVals = (1:vsxParams.P.numRays).*lWidth; 
zVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,size(iq,2));%


figure
ax1 = axes;
B80 = bmode(iq',80);
imagesc(ax1,xVals,zVals,B80);colormap(ax1,'gray');
hold on 

%% Load features 
load([topDir,'/Features',num2str(iImage),'.mat'])

ax2 = axes;
colorData = B80;
imagesc(ax2,xVals,zVals,colorData,'AlphaData',mask)
colormap(ax2,'summer')
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

kLength_BSC_samples = 120;
kLength_LENGTH = 0.5 * kLength_BSC_samples/samplesPwavel*lambda;


[~, zIdx] = min(abs(zVals - zSelect*1e-3));
axIdxs = zIdx-round(kLength_BSC_samples/2) : zIdx + round(kLength_BSC_samples/2) -1 ;
    


%% SUM index ROC testing 


% for all sumIdxs{

%   Pick EML
    
%   for points in EML{
%       Do segmentation (provided segmentation percentage is new) 
%       Compute specificity and sensitivity}
    
%   plot ROC curve 
    
%   compute area under ROC curve}

% plot ROC area as a function of sumIdx 

% determine which maxims the ROC area 

kWidth = 5;
oLap = 0.8;


load([speckleDir,'/COSIEoutput_adaptive/COSIEoutput',num2str(1),'.mat'])


allLines = testData.rayIdxs;
nPossibleKernels = size(idxClustering(allLines,kWidth,oLap),2);


mask2 = mask(axIdxs,allLines);
mask3 = any(mask2,1);

[manSegKernels,manSegKernelIDs]  = idxClustering(find(~mask3),kWidth,oLap);


for iSumIdx = 32

    load([speckleDir,'/COSIEoutput_adaptive/COSIEoutput',num2str(iSumIdx),'.mat'])
    
    segPct2 = 100;
    
    vble = sum(testData.cohAll(:,1:iSumIdx),2);

    specif = zeros(1,size(EML,2));
    sens   = zeros(1,size(EML,2));
    
    for iPoints = 1:size(EML,2)
        
        segBool = vble > EML(1,iPoints) & vble < EML(2,iPoints);
        
        segIdxs = find(segBool);
        
        [kIdxs , kIdxs2] = idxClustering(segIdxs,kWidth,oLap);
       
        nKernels = size(kIdxs,2);
        segPct =  100*(1- nKernels/nPossibleKernels);
        
        %values of kIdxs2 that are in manSegKernelIDs are true positives
        TP  = intersect(kIdxs2,manSegKernelIDs);
        
        %values of kIdxs2 that are not in manSegKernelIDs are false positives
        FP  = setdiff(kIdxs2,manSegKernelIDs);
        
        %values in the set 1:95 that are not in kIdxs2 AND not in
        %manSegKernelIds are true negative 
        TN = intersect(setdiff(1:95,kIdxs2),setdiff(1:95,manSegKernelIDs));

        %values that are in the set manSegKernelIDs but not in kIdxs2 are
        %false negatives 
        FN =  setdiff(manSegKernelIDs,kIdxs2);
       
    end
end





