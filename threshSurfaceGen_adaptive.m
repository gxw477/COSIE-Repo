function [bscSurface,segmentSurface] = threshSurfaceGen_adaptive(cohValues,bscValues,params)
%threshSurfaceGen Generates BSC surface in coherence threshold space
%   cohValues : [1 x nlines] vector of coherence values
%   bscValues : [1 x nlines] vector of bsc values 
%   params.dTH : coherence threshold increment
%   params.APSize : aperture size used
%

    nLines = size(cohValues,2);

    [s1, s2] = size(cohValues);
    [s3, s4] = size(bscValues);

    if s1~= 1 && s2 ~= 1 
        error('coherence VECTOR (1xn) please')
    end
    if s3~= 1 && s4 ~= 1 
        error('bsc VECTOR (1xn) please')
    end
    if s1 > s2 
        cohValues = cohValues';
        [s1, s2] = size(cohValues);
    end
    if s3 > s4 
        bscValues = bscValues';
        [s3, s4] = size(bscValues);

    end
    if s2 ~= s4
        error('Same length vectors please')
    end
    
    nLines = s4;
    
    %Preserves order in cohValues whilst making an ordered vector for
    %analysis
    uTH = sort(cohValues);
    lTH = uTH;
    
    %segBool = zeros(nLines,length(lTH),length(uTH));
    bscSurface = zeros(length(lTH),length(uTH));
    segmentSurface = zeros(length(lTH),length(uTH));


    
     
    for ilTH = 1:length(lTH)
        for iuTH = ilTH:length(uTH)
            
            boolVec = cohValues > lTH(ilTH) & cohValues <= uTH(iuTH);
            %segBool(:,ilTH,iuTH) = boolVec;
            
            sValue = bscValues(boolVec);
            
            if isempty(sValue)
                sValue = 0;
            end
            
            bscSurface(iuTH,ilTH) = abs(mean(sValue));
            segmentSurface(iuTH,ilTH) = length(find(boolVec))/length(boolVec);
 
        end
    end
    
    %bscSurface = bscValues;

end