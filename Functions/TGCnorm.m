

function [imgOut] = TGCnorm(imgIn,vsxStruct,yVals)

    
      
    
    
    %input TGC: defined in manual
    iLine = 0:1023;
    %output gain
    oLine = linspace(0,40,length(iLine));
    
    %gain as defined by TGC waveform
    gainVals = interp1(iLine,oLine,cPtVector);

    gainValsNorm = gainVals-max(gainVals);
    gainValsAdd =  -gainValsNorm;

    dBmultiply = 10.^(gainValsAdd./20);
    
    imgOut = imgIn(:,1:nSamples2).*dBmultiply';

end