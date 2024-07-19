
function [mEdgeSpect ,y0 ] = edgeDetectionMulti(fullIM,fs,yVals,kLengthPER,c0,f0,xBool)

    %get fast time axis as the first dimension
    if size(fullIM,2)>size(fullIM,1)
        fullIM = fullIM';
    end

    figure 
    iq = rf2iq(fullIM(:,xBool),fs);
    bm = bmode(iq,50);
    
    xDummy = 1:size(iq,1);

    imagesc(1:size(iq,2),yVals,bm)
    colormap gray
    ylim([0 10e-3])
    
    kLength_LENGTH = 0.5*kLengthPER*(1/f0)*c0;

    y0 = input('Interface Depth')

    [ ~ , yIdx] = min(abs(yVals-y0));
    
    samplesPwavel = 4;
    kLengthSamples = kLengthPER*samplesPwavel;
    win = hanning(kLengthSamples);

    yBool = (yIdx - kLengthSamples/2):(yIdx+kLengthSamples/2-1);


    edgeRegion = fullIM(yBool,xBool).*win;
    spectVals = abs(fft(edgeRegion));

    df = fs/(kLengthSamples-1);
    fVals = (0:kLengthSamples-1).*df;

    [~ , nF] = min(abs(fVals-f0));

    mEdgeSpect = mean(spectVals(nF,:));
    

end