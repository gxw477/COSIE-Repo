

function [EML ] = emlGen_adaptive(threshSurf,mubsGT,thVector)
% emlGen generates eml based on ground truth speckle value 
%       threshSurf  : threshold surface in coherence space NxN 
%       mubsGT      : ground truth speckle value 1
%       thVector    : defines co-ordinated for threshSurf 1xN
%       ~
%       EML         : Error minimising line (2xN)
    

    %F = griddedInterpolant(thVector,thVector);

    errMubs = abs(threshSurf-mubsGT);

    N = length(thVector);

    EML = zeros(2,N);

    for iTH = 1:N

        [~,errMIN_IDX] = min(errMubs(:,iTH));
        
        
        if errMIN_IDX < iTH
            errMIN_IDX = nan;
            EML(:,iTH) = [nan,nan];
            
        else
            EML(:,iTH) = [thVector(iTH),thVector(errMIN_IDX)];

        end
    end





end