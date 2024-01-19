
clear 
close all 

%Image Analysis

topDir = 'PhantomExperimentsL74_QuadInterp/QA040GSE/pt7dB/';

%saveDir = 'PhantomExperimentsL74_QuadInterp/QApht_angle/';
saveDir = topDir;

vsxParams = load([topDir,'/VSXoutput.mat']);



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

[~, bscFocIdx] = min(abs(rVals - vsxParams.TX(1).focus*lambda));



%% Define coherence kernel properties 

%kernel length (wavels)
nWavels = 3;
%kernel length (samples)
kLength_COH = round(vsxParams.Receive(1).samplesPerWave);



  
%% 
nImages =  12;
kIdxs = cell(nImages,1);

axIdxs = bscFocIdx-round(kLength_BSC_samples/2) : bscFocIdx + round(kLength_BSC_samples/2) -1 ;

fs = vsxParams.Trans.frequency*1e6*vsxParams.Receive(1).samplesPerWave;

fVals = (0:length(axIdxs)-1).*fs/(length(axIdxs));
df = fVals(2)-fVals(1);
nF = round(vsxParams.Trans.frequency*1e6/df);
 

kCount = 1; 
rayCount= 0; 

win = tukeywin(120,0.1)';


for iImage = 1:nImages

    iImage

    load([topDir,'BFimgData',num2str(iImage),'.mat'])
   
    iq2 = iq(:,axIdxs);

    figure
    imagesc(log(abs(fullIM')),[2.5 11]); colormap gray 
    hold on 
    plot([1 size(fullIM,1)],bscFocIdx.*[1 1],'-','color','red')
    plot(120.*[1 1], [axIdxs(1) axIdxs(end)],'-','color','red')
    
    
    colormap gray

    pause(0.5)

    lRl = 17;%input('Left ray line: ');
    rRl = 111;%input('Right ray line: ');

    rayIdxs = lRl:rRl;

    if length(axIdxs) ~= kLength_BSC_samples
        error('check bsc kernel length')
    end

    %% BSC + Coherence Calc

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

