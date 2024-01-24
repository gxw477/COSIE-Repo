
function bscEstimate = COVsegmentation(vble,EML,spectAll,kWidth,oLap,nPossibleKernels)
    
    segPctIdx = 2;
    segPct = 0;
    
    kIdxs = idxClustering(1:size(spectAll,1),kWidth,oLap);
    nKernels = size(kIdxs,2);

    bscValuesIK = zeros(1,nKernels);

    for iKernels = 1:nKernels
        bscValuesIK(iKernels) = abs(mean(spectAll(kIdxs{iKernels})));
    end

    bscEstimate(:,1) = [mean(bscValuesIK) ;std(bscValuesIK); 0; skewness(bscValuesIK);kurtosis(bscValuesIK)];

    
    while segPct < 50
    
        segBool = vble > EML(1,segPctIdx) & vble < EML(2,segPctIdx);
        
        segIdxs = find(segBool);
        
        kIdxs = idxClustering(segIdxs,kWidth,oLap);
    
        nKernels = size(kIdxs,2);
    
        bscValuesIK = zeros(1,nKernels);
    
        for iKernels = 1:nKernels
            bscValuesIK(iKernels) = abs(mean(spectAll(kIdxs{iKernels})));
        end
    
        segPct =  100*(1- nKernels/nPossibleKernels);
    
        bscEstimate(:,segPctIdx) = [mean(bscValuesIK) ;std(bscValuesIK);segPct; skewness(bscValuesIK);kurtosis(bscValuesIK)];
        
        segPctIdx = segPctIdx+1;
    
    end    
    
end
