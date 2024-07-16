

function [ channels] = getElementData(M,SIG,methodIdx)

%getElementData Takes beamforming matrix and signal, and outputs to give
%channel specific data for coherence calculations. Current Version 
%works for any interpolation method
%   Input
%   ~~~~~~~~~~~~~~~~~~
%   M: Beamforming matrix 
%   SIG: Signal
%   methodIdx: method of interpolation
%
%      1: 'nearest'   - nearest neighbor interpolation
%      2: 'linear'    - (default) linear interpolation
%      3: 'quadratic' - quadratic interpolation
%      4: 'lanczos3'  - 3-lobe Lanczos interpolation
%      5: '5points'   - 5-point least-squares parabolic interpolation
%      6: 'lanczos5'  - 5-lobe Lanczos interpolation
%
%
%   Output
%   ~~~~~~~~~~~~~~~~~~
%   bfElement : Beamformed Element data
%
%

    if nargin == 0
        %load('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\MUST_GEORGE\sparseTestData.mat')
    end
    
    %   SIG has size  [nChannels*(nSamples+1) x 1M" =
    %   M has size [nSamples x length(SIG)]
    %   beamformed line would be bFSIG = M*SIG
    if size(SIG,2) > 1
        error('You have made a mistake my friend, check SIG size on passing to getElementData, please and thankyou')
    end
    
    nSamples = size(M,1);
    nChannels = length(SIG)/(size(M,1));
    
    channels = zeros(nSamples,nChannels);

    if nChannels > 128 
        nChannels
    end

    lVec = zeros(1,nSamples);

    for iSample = 1:nSamples
        
        %add small amount to signal to avoid zero values being ignored in
        %next line 
        SIG = SIG + 1e-11; 
        
        %ith pixel calc. Selects data samples for summation and weights to
        %do the linear interp
        vec = nonzeros(M(iSample,:).*SIG');
        vec2 = vec(1:methodIdx:end);

        %sum the scaled values
        for i = 2:methodIdx
            size(vec(1:methodIdx:end));
            size(vec(i:methodIdx:end));

            vec2 = vec2 + vec(i:methodIdx:end); 

        end
        
        %subtract off small amount (n samples so multiplied by n).
        vec2 = vec2 - 1e-11*methodIdx;
        
        channels(iSample,1:length(vec2)) = vec2; 

        lVec(iSample) = length(vec);


        if lVec(iSample) < 128
            vec = [vec;zeros(128-lVec(iSample),1)];
        end

        if 0 %lVec(iSample) > 128
            iSample
        end
    end

end


