clear 

close all 

topDir = [uigetdir,'\'];
cd(topDir)

%fileIdxs = 5:5:50;
fileNames = ls('Liver*');
%fileNames = fileNames(3:end,:);

%load(uigetfile)

SOS = input('Speed of sound :')

for iFrame = 1:size(fileNames,1)

    VSXfileOption = 2 ;
    
    if VSXfileOption == 1
        
        load([topDir,'VSXinit.mat'])
        samplesPerAcq = (Receive(1).endDepth-Receive(1).startDepth)*4;
    
    elseif VSXfileOption == 2
        load([topDir,'VSXoutput.mat'])
        samplesPerAcq = Receive(1).endSample - Receive(1).startSample + 1;
         
    end
    
    
    %load([topDir,'vsxResult_newSample.mat'])
    %load([topDir,'Img_',num2str(fileIdxs(iFrame)),'.mat'])
    %load(fileNames(iFrame,:))
    %close all

    load([fileNames(iFrame,:)])
    
    RfFrame1 = double(RcvData{1}); 
    
    
    clearvars RfRaw
    
    %% Calculate reconstruction positions
    
    %lambda = (Resource.Parameters.speedOfSound/(Trans.frequency*1e6));
    lambda = SOS/(Trans.frequency*1e6);

    samplesPwavel = Receive.samplesPerWave; 
     
    fs = Trans.frequency*1e6*samplesPwavel;
    
    dY_LAMBDA = 1/samplesPwavel/2;
    dY_M = lambda/samplesPwavel/2;
    
    %elPosX = Trans.ElementPos(:,1);
    %elPosZ = Trans.ElementPos(:,3);
    
    
    bformX = zeros(1,samplesPerAcq); 
    
    samplesPwavel = Receive(1).samplesPerWave; 
    
    
    %% Varargin 
    
    %params.fs = fs;
    elPos1 = [ Trans.ElementPos(1,1) Trans.ElementPos(1,2) Trans.ElementPos(1,3)];
    elPos2 = [ Trans.ElementPos(2,1) Trans.ElementPos(2,2) Trans.ElementPos(2,3)];
    
    if Trans.units == 'wavelengths'
        factor = lambda;
    else 
        factor = 1e-3;
    end
    
    params.pitch = Trans.spacingMm*1e-3;
    params.fc = Trans.frequency*1e6;
    params.kerf = [];
    %params.width = Trans.elementWidth *factor;
    params.fs = params.fc *samplesPwavel;
    params.Nelements = Trans.numelements;
    %params.bandwidth = 79;
    %params.radius = Trans.radiusMm*1e-3;
    params.height = Trans.elevationApertureMm*1e-3;
    params.focus = Trans.elevationFocusMm*1e-3;
    params.c = SOS; %Resource.Parameters.speedOfSound;
    params.t0 = 0; %Receive(1).startDepth*lambda*2/params.c  ;

    %% compute DAS for each line
    
    fullIM = zeros(P.numRays, samplesPerAcq);
    
    %period of f0
    tau = 1/(Trans.frequency*1e6);
    
    %Pre-allocate memory for channel stack variable 
    
    apLengths = zeros(1,P.numRays);
    
    for i = 1:P.numRays
        
        apLengths(i) = length(find(TX(i).Apod));
    
    end 
    
    RfFrame2 = restructureFrame(RfFrame1,Trans.Connector);
    
    clearvars RfFrame1
    
    channelCheck = sum(RfFrame2,1);
    omitChannel = channelCheck == 0;
    
    channelStack = zeros(P.numRays, samplesPerAcq , max(apLengths));
    
    if VSXfileOption ==1 
        iPix = 1; 
    end
    
    
    for i = 1:P.numRays
        
        i
        
        [~, cIdx ] = max(TX(i).Delay);
    
        tDelay = (TX(i).Delay) .*tau;
        %omitBool = tDelay == 0;
        omitBool = TX(i).Apod == 0;
      
        tDelay(omitBool) = nan;
    
        %apElPos = [Trans.ElementPos(~omitBool,1), Trans.ElementPos(~omitBool,3)];
        %chordLength = sqrt((apElPos(1,1)-apElPos(end,1))^2 +(apElPos(1,2)-apElPos(end,2))^2);
        %chordAngle = 2*asin(chordLength/(2*Trans.radius));
        %dz = Trans.radius * (1-cos(chordAngle/2))*lambda;
            
        acqI = RfFrame2(:,~omitChannel);
        acqI = acqI(Receive(i).startSample:Receive(i).endSample,~omitBool);
        
        samplesPerAcq2 = Receive(i).endSample - Receive(i).startSample + 1;
        
        if samplesPerAcq2 ~= samplesPerAcq
            error('Samples per acqusition error')
        end
        
        bformY =  linspace(Receive(i).startDepth,Receive(i).endDepth,samplesPerAcq2);
    
        params.Nelements  = length(tDelay(~omitBool));
        [bfSig,channelData,~] = mustUGeorge(acqI, zeros(1,length(bformY)), bformY.*lambda,tDelay(~omitBool),params,'quadratic');
         
        fullIM(i,:) = bfSig;
       
        channelStack(i,:,1:size(channelData,2)) = channelData;
        
    end
    
    
    iq = rf2iq(fullIM,params.fs);
    
    I = IData{1};
    Q = QData{1};
    IQ = abs(I + 1i.*Q);
    
    
    %% Calculate Coherence
    
    
    %kernel length (wavels)
    nWavels = 3;
    %kernel length (samples)
    kLength = 120;%round(samplesPwavel*nWavels); 
    
    kLength_LENGTH = (bformY(1+kLength)-bformY(1))*lambda;
    nKernels = floor(samplesPerAcq/kLength);
    
    kY = min(bformY)*lambda + kLength_LENGTH/2 + (0:nKernels-1).*kLength_LENGTH;
    
    RMat = zeros(P.numRays,nKernels,max(apLengths));
    
    tic 
    for iRay = 1:P.numRays
    
    iPix = 1;
        
        for iK = 1:nKernels
            
            scatData = squeeze(channelStack(iRay,iPix:(iPix + kLength-1),:));
            RMat(iRay,iK,:) = CoherenceAnalysisFN(scatData);
            iPix = iPix + kLength;
            
        end
    
    end
    toc
    
    cohValues = sum(RMat,3);
    
    %% Multiply coherence by B-mode
    
    [xQbm ,yQbm ] = ndgrid(1:P.numRays,bformY.*lambda);
    [xQcoh,yQcoh] = ndgrid(1:P.numRays,kY);
    
    F = griddedInterpolant(xQcoh,yQcoh,cohValues);
    
    CoherenceQ = F(xQbm,yQbm);
    
    %Coherence weighted IQ image
    compImage = CoherenceQ.*iq;
    
    %% Plots 
    
    yVals = lambda.*linspace(Receive(1).startDepth,Receive(1).endDepth,samplesPerAcq);%
    
    %[rMeshed, thMeshed] = meshgrid(rVals,Angle);
    %[zPolar , xPolar  ] = pol2cclear 
    % art(thMeshed,rMeshed);
    
    startIdx = 250;
    endIdx = 1500; 
    %zPolar = zPolar - radius*lambda; 
    %zPolar2 = zPolar(:,startIdx:endIdx);
    %xPolar2 = xPolar(:,startIdx:endIdx);.


    B = bmode(iq(:,startIdx:endIdx),20);
        
    save([topDir,'\BFimgData',num2str(SOS),'.mat'])

    clearvars -except iFrame topDir fileNames

end
