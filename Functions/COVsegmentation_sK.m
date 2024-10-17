
function bscEstimate = COVsegmentation_sK(vble,EML,spectAll,kWidth,oLap)
    
    segPctIdx = 1;
    segPctL = 0;
    
    nPossibleKernels = size(idxClustering(1:length(vble),kWidth,oLap),2);    

    kIdxs = idxClustering(1:length(spectAll),kWidth,oLap);
    nKernels = size(kIdxs,2);

    bscValuesIK = zeros(1,nKernels);

    for iKernels = 1:nKernels
        bscValuesIK(iKernels) = mean(abs(spectAll(kIdxs{iKernels})));
    end

    %Unsegmented Values
    bscEstimate(:,1) = [mean(bscValuesIK) ;std(bscValuesIK); 0; 0; skewness(bscValuesIK);kurtosis(bscValuesIK)];

    
    for segPctIdx = 1:length(EML(1,:))
    
        segBool =   vble > EML(1,segPctIdx) & vble < EML(2,segPctIdx);
       
        
        segIdxs = find(segBool);
        
        kIdxs = idxClustering(segIdxs,kWidth,oLap);
    
        nKernels = size(kIdxs,2);
    
        bscValuesIK = zeros(1,nKernels);
    
        for iKernels = 1:nKernels
            bscValuesIK(iKernels) = abs(mean(spectAll(kIdxs{iKernels})));
        end
    
        segPctL =  100*(1-length(find(segBool))/length(segBool));
        segPctK =  100*(1- nKernels/nPossibleKernels);
    
        bscEstimate(:,segPctIdx+1) = [mean(bscValuesIK) ;std(bscValuesIK);segPctL; segPctK; skewness(bscValuesIK);kurtosis(bscValuesIK)];
        
        %segPctIdx = segPctIdx+1;
    
    end    
    
end
