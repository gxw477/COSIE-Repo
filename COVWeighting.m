function [bscEstimate] = COVWeighting(inVble,binCentre,pDist,powf0,kWidth,oLap)

%   COVWeighting : Weights a ray line based on its posiiton within the
%   probability distribution of the variable for the speckle ground truth data
%   set
%   
%   inVble : coherence or texture 
%   pDist  : probability distribution of the variable
%
%   bscEstimate : 5x2 vector : Mean , SD , seg %, skewness, kurtosis


    bscEstimate = zeros(5,2);
        
    pPoints = interp1(binCentre,pDist,inVble);
    pPoints(isnan(pPoints))= 0;

    
    %groups ray lines based on kernel width and overlap
 
    nPossibleKernels = size(idxClustering(1:length(inVble),kWidth,oLap),2);

    kIdxs1 = idxClustering(1:length(powf0),kWidth,oLap);
    nKernels1 = size(kIdxs1,2);
    bscValuesIK = zeros(1,nKernels1);


    for iKernels = 1:nKernels1
        bscValuesIK(iKernels) = abs(mean(powf0(kIdxs1{iKernels})));
    end

    bscEstimate(:,1) = [mean(bscValuesIK) ;std(bscValuesIK); 0; skewness(bscValuesIK);kurtosis(bscValuesIK)];


    idxsAll = 1:length(inVble);
    idxsAll = idxsAll(pPoints>0);

    kIdxs2 = idxClustering(idxsAll,kWidth,oLap);
    
    %number of kernels total 
    nKernels = size(kIdxs2,2);
    bscValuesIK = zeros(1,nKernels);

    for iKernel = 1:nKernels
        
        normFactor = 1/(sum(pPoints(kIdxs2{iKernel})));
        weightFactor = normFactor.*pPoints(kIdxs2{iKernel});

        if size(weightFactor,2)>size(weightFactor,1)
            weightFactor = weightFactor';
        end

        %bsc value per kernel
        bscValuesIK(iKernel) = abs(mean(weightFactor.*powf0(kIdxs2{iKernel})));
    end


    segPct =  100*(1- nKernels/nPossibleKernels);
    
 
    %mean value across kernels
    bscEstimate(:,2) = [mean(bscValuesIK);std(bscValuesIK);segPct;skewness(bscValuesIK);kurtosis(bscValuesIK)];   

end