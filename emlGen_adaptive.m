

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

    %Tolerance1 : How close can the upper and lower thresholds be ? (units
    %of indices of thVector
    tol1 = 10;

    %Tolerance 2 : How much can the EML jump by (max) ? 
    tol2 = 1;


    for iTH = 1:N

        [~,errMIN_IDX] = min(errMubs(:,iTH));
        
        
        if (errMIN_IDX - iTH) < tol1 && iTH == 1

           
            errMIN_IDX = nan;
            EML(:,iTH) = [nan,nan];
        elseif  (errMIN_IDX - iTH) < tol1 &&  abs(thVector(errMIN_IDX)- thVector(errMIN_IDX)) > tol2 
            errMIN_IDX = nan;
            EML(:,iTH) = [nan,nan];
        else
            EML(:,iTH) = [thVector(iTH),thVector(errMIN_IDX)];

        end
    end





end