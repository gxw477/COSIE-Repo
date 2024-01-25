function [spectOut] = spectAveraging(spectIn, kWidth, oLap)
%[spectOut,] = spectAveraging(spectIn, kWidth, oLap)
%           spectIN :   Original 
%           kWidth  :   kernelWidth(lines)
%           oLap    :   overlap betweeen kernels (decimal)
%
%           spectOut:   Kernel averages
%           cents   :   kernel centres

    
    nLines = size(spectIn,1);
    nSamples = size(spectIn,2);

    nKernels = floor(nLines/kWidth) + 1 ;

    %flap = nLines - kWidth*nKernels

    shift = round((1-oLap)*kWidth); 
    
     
    spectOut = zeros(nKernels,nSamples);
    
    for iK = 1:nKernels 
        
        if iK < nKernels
            left = (iK-1)*shift + 1;
            right = left+ kWidth -1; 
        else
            right = size(spectIn,1);
            left = right - kWidth;
        end

        spectOut(iK,:) = mean(spectIn(left:right,:),1);
        
    end

end
