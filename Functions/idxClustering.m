function [kIdxs, kIdxs2] = idxClustering(idxs,width,oLap)
%idxClustering finds idxs within # of raylines `width' and overlap `oLap'

    kIdxs = cell(1,1);
    
    kN = 1;

    iRayline = 3;
    sepn = round(oLap*width);
    
    endRayline = width;
    %startRayline = iRayline - 2;

    while endRayline <= max(idxs) 
        
        vec = (endRayline-width + 1):endRayline;

        compBool = ismember(vec,idxs);
        
        if all(compBool)
            kIdxs{kN} = vec;
            kN = kN+1;
        end
        
        endRayline = endRayline + width - sepn;

    end
    
    kIdxs2 = unique(cell2mat(kIdxs));

end