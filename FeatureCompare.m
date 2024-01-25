


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

depthSelect = input('Depth of interest (mm) : ');

speckleDir = [speckleDir,'Z',num2str(depthSelect),'\'];
topDir2 = [topDir,'Z',num2str(depthSelect),'\'];
load([topDir2,'COSIEinput',num2str(iImage),'.mat'])


tgcBool = isequal(vsxParams.TGC,vsxParams2.TGC)

if ~tgcBool
    error('TGC altered between Speckle and test image')
end

figure
B80 = bmode(iq',80);
imagesc(B80);colormap gray

%% Load features 
load([topDir,'/Features',num2str(iImage),'.mat'])


hold on 
for i = 1:size(xFeature2,1)
    plot([xFeature2(i,1) xFeature2(i,2)],yFeature2(i).*[1 1],'-','Color','red')
    hold on
end


%%





