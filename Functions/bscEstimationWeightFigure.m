function [] = bscEstimationWeightFigure(bscSpeckleBf,bscSpeckleSTD,bscEstimate_ENV_WEIGHT,bscEstimate_COH_WEIGHT)

% bscEstimationWeightFigure(bscSpeckleBf,bscSpeckleSTD,bscEstimate_ENV_WEIGHT,bscEstimate_COH_WEIGHT)

% bscSpeckleBf : Ground truth speckle 
% bscSPeckleSTD : Ground truth speckle uncertainty
% bscEstimate_ENV_WEIGHT : array containing weight estimates and seg %
% using SNR
% bscEstimate_COH_WEIGHT : array containing weight estimates and seg %
% using Coherence
%

figure
errorbar(-10,bscSpeckleBf,bscSpeckleSTD,'r.')
hold on 
errorbar(0,bscEstimate_ENV_WEIGHT(1,1),bscEstimate_ENV_WEIGHT(2,1),'k.')

errorbar(10,bscEstimate_COH_WEIGHT(2,1),bscEstimate_COH_WEIGHT(2,2),'k.')
errorbar(20,bscEstimate_COH_WEIGHT(2,1),bscEstimate_COH_WEIGHT(2,2),'k.')

xlim([-15 30])
xticks(-10:10:20)
xticklabels({'Ground Truth','Unseg','Coh weights','SNR weights'})

title('Weighting Results')



end