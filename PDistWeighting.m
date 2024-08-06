function [bscEstimate] = PDistWeighting(inVble,binCentre,pDist,powf0,kWidth,oLap)

%   PDistWeighting : Weights a ray line based on its posiiton within the
%   probability distribution of the variable for the speckle ground truth data
%   set
%   
%   inVble : coherence or texture 
%   pDist  : probability distribution of the variable
%
%   bscEstimate : 5x2 vector : Mean , SD , seg %, skewness, kurtosis


    bscEstimate = zeros(6,2);
        
    weights = interp1(binCentre,pDist,inVble);
    
    %weights = ones(1,length(inVble));%interp1(binCentre,pDist,inVble);
    
    weights(isnan(weights))= 0;
    
    %sWeights = sum(weights);
    
    %groups ray lines based on kernel width and overlap
 
    nPossibleKernels = size(idxClustering(1:length(inVble),kWidth,oLap),2);

    kIdxs1 = idxClustering(1:length(powf0),kWidth,oLap);
    nKernels = size(kIdxs1,2);
    bscValuesIK = zeros(1,nKernels);


    for iKernels = 1:nKernels
        bscValuesIK(iKernels) = abs(mean(powf0(kIdxs1{iKernels})));
    end

    bscEstimate(:,1) = [mean(bscValuesIK) ;std(bscValuesIK); 0 ; 0; skewness(bscValuesIK);kurtosis(bscValuesIK)];

    bscValuesIK = zeros(1,nKernels);
    kWeightsIK = ones(1,nKernels);

    for iKernel = 1:nKernels
        
        pIK = weights(kIdxs1{iKernel});
        spIK = sum(pIK);
        kWeightsIK(iKernel) = spIK;
   
        %bsc value per kernel
        bscValuesIK(iKernel) = abs(mean(powf0(kIdxs1{iKernel})));

       
    end

    %Adds up the weigthing instead of the eliminated kernel amount
    segPctL =  100*(1-length(find(weights))/length(inVble));       
    segPctK =  100*(1- nKernels/nPossibleKernels);
    
    spIK2 = sum(kWeightsIK);

    %mean value across kernels
    bscEstimate(:,2) = [sum(kWeightsIK.*bscValuesIK)./spIK2;std(bscValuesIK,kWeightsIK);segPctL;segPctK;skewness(bscValuesIK);kurtosis(bscValuesIK)];   

end