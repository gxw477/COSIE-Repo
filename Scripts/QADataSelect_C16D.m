
clear 
close all 

%Image Analysis

%topDirMaster = [uigetdir('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\','Select Analysis Folder'),'\'];
%topDirMaster = ['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\G218L74_2\']%[uigetdir('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74\','Select Analysis directory'),'\'];
topDirMaster = [uigetdir,'\'];

vsxParams = load([topDirMaster,'\VSXoutput.mat']);



%% Define BSC kernel properties
lambda = vsxParams.lambda;

%dtheta = vsxParams.Angle(2)-vsxParams.Angle(1);

kWidth_BSC_lines = 5;
kLength_BSC_samples = 120;
oLap = .80;

%lWidthF = (vsxParams.Trans.radius + vsxParams.TX(1).focus)*lambda*dtheta;

samplesPerAcq = vsxParams.Receive(1).endSample - vsxParams.Receive(1).startSample + 1;
   

%already defined, but to remind you
rVals = lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,samplesPerAcq)  ;
lWidth = (vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1))*lambda;

xVals = vsxParams.Angle;


[~, bscFocIdx] = min(abs(rVals - vsxParams.TX(1).focus*lambda));



%% Define coherence kernel properties 

%kernel length (wavels)
nWavels = 3;
%kernel length (samples)
kLength_COH = round(vsxParams.Receive(1).samplesPerWave);



  
%% 
fNames = ls([topDirMaster,'\BFimgDataTGCcorr\BFImg*']);
nImages = size(fNames ,1);


kIdxs = cell(nImages,1);


fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;

fVals = (0:kLength_BSC_samples-1).*fs/(length(kLength_BSC_samples));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);

kCount = 1; 
rayCount= 0; 

wOption = input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');

if wOption == 1 
    win = ones(1,kLength_BSC_samples);
    wName = 'Rect';
elseif wOption == 2
    win = tukeywin(kLength_BSC_samples,0.1)';
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welch';
end

for i = 1:128 
    n(i) = sum(vsxParams.TX(i).Apod);
end

rayIdxs = find(n == max(n));
%sumIdx = max(n)

sumImg = zeros(nImages,2048);

for iImage = 1:nImages

    iImage
    
    load([topDirMaster,'\BFimgDataTGCcorr\',fNames(iImage,:)])
    
    figure
    imagesc(xVals,yVals,log(abs(fullIM')),[2.5 11]); colormap gray 
    hold on 
    plot([xVals(1) xVals(end)],vsxParams.TX(1).focus*lambda.*[1 1],'-','color','red')
    colormap gray
    
    fullIM2 = TGCnorm(fullIM,vsxParams,yVals);

    sumImg(iImage,:) = sum(log(abs(fullIM)),1);

end



