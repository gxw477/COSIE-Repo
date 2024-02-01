function [kIdxs, kIdxs2] = idxClustering(idxs,width,oLap)
%idxClustering finds idxs within # of raylines `width' and overlap `oLap'

    kIdxs = cell(1,1);
    
    kN = 1;

    iRayline = 3;
    sepn = oLap*width;
    
    if sepn ~= round(sepn)
        error('Check overlap and kernel width')
    end

    endRayline = iRayline + 2;
    %startRayline = iRayline - 2;

    while endRayline <= max(idxs) 
        
        vec = (endRayline-width + 1):endRayline;

        compBool = ismember(vec,idxs);
        
        if all(compBool)
            kIdxs{kN} = vec;
            kIdxs2(kN:(kN + length(vec) -1)) = vec;
            kN = kN+1;
        end
        
        endRayline = endRayline + width - sepn;

    end

end