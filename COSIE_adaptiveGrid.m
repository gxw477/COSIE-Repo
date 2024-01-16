
function [bscSurface,segSurface,EML,pctSeg1,redEML,pctSeg2,params] = COSIE_adaptiveGrid(cohValues,bscValues,params)
    % [bscSurface,segSurface,EML,pctSeg1,redEML,pctSeg2,params] = COSIE(cohValues,bscValues,params)
    %   cohValues : [1 x nlines] vector of coherence values
    %   bscValues : [1 x nlines] vector of bsc values 
    %   params.APSize : aperture size used
    

    [bscSurface,segSurface] = threshSurfaceGen_adaptive(cohValues,bscValues,params);
    
    %order not needed in preserving indices for the EML (I think).
    thVector = sort(cohValues);

    EML  = emlGen_adaptive(bscSurface,abs(mean(bscValues)),thVector);
   
    [ redEML,pctSeg1,pctSeg2 ] =  emlClassify(cohValues,EML);
    
end