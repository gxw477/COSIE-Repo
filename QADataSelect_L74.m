
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



for iImage = 1:7

    iImage

    load([topDirMaster,'BFimgData',num2str(iImage),'.mat'])
   
    
    figure
    imagesc(xVals,yVals,log(abs(fullIM')),[2.5 11]); colormap gray 
    hold on 
    plot([xVals(1) xVals(end)],vsxParams.TX(1).focus*lambda.*[1 1],'-','color','red')
    colormap gray

    lRl = 17;
    rRl = 111;
    rayIdxs = lRl:rRl;

    pause(0.5)

    for zVals = 15:5:40

        zSelect = zVals*1e-3;
        [~, zIdx] = min(abs(rVals - zSelect));
        axIdxs = zIdx-round(kLength_BSC_samples/2) : zIdx + round(kLength_BSC_samples/2) -1 ;
        
        if length(axIdxs) ~= kLength_BSC_samples
            error('check bsc kernel length')
        end
    
        plot([xVals(1) xVals(end)],zSelect.*[1 1],'-','color','green')
        plot(120*lWidth.*[1 1],[rVals(axIdxs(1)) rVals(axIdxs(end))],'-','color','green')
        
        saveDir = [topDirMaster,wName,'/Z',num2str(round(zSelect*1e3)),'/'];

        if ~exist(saveDir)
            mkdir(saveDir)
        end
        
        bscLines = fullIM(rayIdxs,axIdxs);
    
        if wOption == 1 || wOption == 2
            winMatrix = ones(size(bscLines)).*win;
            bscLines = fullIM(rayIdxs,axIdxs).*winMatrix;
            spect = abs((fft(bscLines,[],2)).^2);
        
        elseif wOption == 3
            h = spectrum.welch;                  % Create a Welch spectral estimator.
            welchObj = psd(h,bscLines','Fs',fs);
            spect = (welchObj.Data').^2; % transpose because spectAveraging takes the vector the other way
        end
    
       
        cohAll = zeros(length(rayIdxs),size(channelStack,3));
        spectAll = zeros(length(rayIdxs),size(spect,2));
        
        for iLine2 = 1:length(rayIdxs)
            
            cohAll(iLine2,:) = CoherenceAnalysisFN(squeeze(channelStack(rayIdxs(iLine2),axIdxs,:)));
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

end

