function [redEML ,pctSeg,pctSegOUT ] = emlClassify(cohValues,EML)
    
    pctSeg= zeros(1,size(EML,2));
    nValues = length(cohValues);
    
    for iTh = 1:size(EML,2)
        
        segBool = cohValues < EML(1,iTh) | cohValues >= EML(2,iTh);
        nSeg = length(find(segBool));
        pctSeg(iTh) = nSeg/nValues*100;
            
    end

    targetSeg = 1:70;  
    redEML = zeros(2,length(targetSeg));
    pctSegOUT = zeros(1,length(targetSeg));

    for iSeg = 1:length(targetSeg)
        
        [~, idx ] = min(abs(pctSeg-targetSeg(iSeg)));
        redEML(:,iSeg) = EML(:,idx);
        pctSegOUT(iSeg) = pctSeg(idx);
        
    end
end