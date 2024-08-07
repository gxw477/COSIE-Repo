
function bscEstimate = COVsegmentation(vble,EML,spectAll,kWidth,oLap,rejectThresh)
    
    segPctIdx = 1;
    segPctL = 0;
    
    nPossibleKernels = size(idxClustering(1:length(vble),kWidth,oLap),2);    

    kIdxs = idxClustering(1:size(spectAll,1),kWidth,oLap);
    nKernels = size(kIdxs,2);

    bscValuesIK = zeros(1,nKernels);

    for iKernels = 1:nKernels
        bscValuesIK(iKernels) = abs(mean(spectAll(kIdxs{iKernels})));
    end

    %Unsegmented Values
    bscEstimate(:,1) = [mean(bscValuesIK) ;std(bscValuesIK); 0; 0; skewness(bscValuesIK);kurtosis(bscValuesIK)];

    
    while segPctIdx < 200
    
        segBool =   vble > EML(1,1) & vble < EML(2,1);
        segBool2 = sum(segBool,2);
        segBool3 = segBool2 > rejectThresh;

        
        segIdxs = find(segBool3);
        
        kIdxs = idxClustering(segIdxs,kWidth,oLap);
    
        nKernels = size(kIdxs,2);
    
        bscValuesIK = zeros(1,nKernels);
    
        for iKernels = 1:nKernels
            bscValuesIK(iKernels) = abs(mean(spectAll(kIdxs{iKernels})));
        end
    
        segPctL =  100*(1-length(find(segBool3))/length(segBool3));
        segPctK =  100*(1- nKernels/nPossibleKernels);
    
        bscEstimate(:,segPctIdx+1) = [mean(bscValuesIK) ;std(bscValuesIK);segPctL; segPctK; skewness(bscValuesIK);kurtosis(bscValuesIK)];
        
        segPctIdx = segPctIdx+1;
    
    end    
    
end
