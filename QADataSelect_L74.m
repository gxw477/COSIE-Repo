
clear 
close all 

%Image Analysis

topDir = [uigetdir('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74_QuadInterp\','Select Analysis Folder'),'\'];


vsxParams = load([topDir,'\VSXoutput.mat']);



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

xVals = (1:128).*lWidth;


[~, bscFocIdx] = min(abs(rVals - vsxParams.TX(1).focus*lambda));



%% Define coherence kernel properties 

%kernel length (wavels)
nWavels = 3;
%kernel length (samples)
kLength_COH = round(vsxParams.Receive(1).samplesPerWave);



  
%% 
nImages =  12;
kIdxs = cell(nImages,1);


fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;

fVals = (0:kLength_BSC_samples-1).*fs/(length(kLength_BSC_samples));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);
 

kCount = 1; 
rayCount= 0; 

win = tukeywin(120,0.1)';


for iImage = 1:nImages

    iImage

    load([topDir,'BFimgData',num2str(iImage),'.mat'])
   
    
    figure
    imagesc(xVals,yVals,log(abs(fullIM')),[2.5 11]); colormap gray 
    hold on 
    plot([xVals(1) xVals(end)],vsxParams.TX(1).focus*lambda.*[1 1],'-','color','red')
    colormap gray

    lRl = 17;
    rRl = 111;
    rayIdxs = lRl:rRl;

    pause(0.5)

    zSelect = input('Select depth of interest: ');
    [~, zIdx] = min(abs(rVals - zSelect));
    axIdxs = zIdx-round(kLength_BSC_samples/2) : zIdx + round(kLength_BSC_samples/2) -1 ;
    
    if length(axIdxs) ~= kLength_BSC_samples
        error('check bsc kernel length')
    end

    plot([xVals(1) xVals(end)],zSelect.*[1 1],'-','color','green')
    plot(120*lWidth.*[1 1],[rVals(axIdxs(1)) rVals(axIdxs(end))],'-','color','green')
    
    saveDir = [topDir,'/Z',num2str(round(zSelect*1e3))];

    %% BSC + Coherence Calc

    if ~exist(saveDir)
        mkdir(saveDir)
    end

    bscLines = fullIM(rayIdxs,axIdxs).*win;

    spect = fft(bscLines,[],2);
    
    cohAll = zeros(length(rayIdxs),size(channelStack,3));
    spectAll = zeros(length(rayIdxs),size(spect,2));
    
    for iLine2 = 1:length(rayIdxs)
        
        cohAll(iLine2,:) = CoherenceAnalysisFN_new(squeeze(channelStack(rayIdxs(iLine2),axIdxs,:)));
        spectAll(iLine2,:) = spect(iLine2,:);

    end
    
    %kCount = kCount + nKsInAverage; 
    %rayCount = rayCount + length(kIdxs{iImage});

    save([saveDir,'\COSIEinput',num2str(iImage),'.mat'],'cohAll','spectAll','fVals','rayIdxs')


    %% Envelope statistics calculations 

    %fullIM
    fullEnv = abs(envelope(fullIM(rayIdxs,axIdxs)));
    
    envMean = mean(fullEnv,2);
    envStd = std(fullEnv,0,2);
    
    save([saveDir,'\EnvStats',num2str(iImage),'.mat'],'envMean','envStd')

    mean(envMean./envStd)

end

