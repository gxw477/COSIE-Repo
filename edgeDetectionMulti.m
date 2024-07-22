
function [mEdgeSpect ,y0] = edgeDetectionMulti(dir)
    
    %kernel length in periods of T0 = 1/f0
    kLengthPER = 10;
    xBool = 17:111;
   
   
    IDF = load([dir,'\IDF_FiltAvg.mat']);

    fnames = ls([dir,'\BFimgData*']);
    nImages = size(fnames,1);

    mEdgeSpect = zeros(1,nImages);
    y0 = zeros(1,nImages);
    
    for iImage = 1:nImages
        
        bfImgData = load([dir,'\BFimgData',num2str(iImage),'.mat']);

        fs = bfImgData.fs;
        c0 = bfImgData.Resource.Parameters.speedOfSound;
        yVals = bfImgData.yVals;
        f0 = bfImgData.Trans.frequency*1e6;
        fullIM = bfImgData.fullIM;
        
        lambda = c0/f0;
        offset = bfImgData.Receive(1).startDepth*lambda;
        wSize = IDF.Filt.wsize*lambda; 

        IDF.Filt.yVals = offset + wSize/2 + (IDF.Filt.wpos )./2*lambda/bfImgData.Receive(1).samplesPerWave;
        
        [~,nF0_IDF] = min(abs(f0/1e6-IDF.Filt.f));

        %get fast time axis as the first dimension
        if size(fullIM,2)>size(fullIM,1)
            fullIM = fullIM';
        end
        
        if 0 
            figure 
            iq = rf2iq(fullIM(:,xBool),fs);
            bm = bmode(iq,50);
            
            xDummy = 1:size(iq,1);
    
            imagesc(1:size(iq,2),yVals,bm)
            colormap gray
            ylim([0 40e-3])
        end

        kLength_LENGTH = 0.5*kLengthPER*(1/f0)*c0;
    
        [~,yIdx]= max(sum(fullIM,2));

        y0_Val = yVals(yIdx);
        [~,y0IDF_idx] = min(abs(y0_Val-IDF.Filt.yVals));

       

        y0(iImage)= y0_Val;
     
        samplesPwavel = 4;
        kLengthSamples = kLengthPER*samplesPwavel;
        win = hanning(kLengthSamples);
    
        yBool = (yIdx - kLengthSamples/2):(yIdx+kLengthSamples/2-1);
        
        edgeRegion = fullIM(yBool,xBool).*win;
        spectVals = abs(fft(edgeRegion));
    
        df = fs/(kLengthSamples-1);
        fVals = (0:kLengthSamples-1).*df;
    
        [~ , nF] = min(abs(fVals-f0));
        
      
        mEdgeSpect(iImage) = mean(spectVals(nF,:))/IDF.Filt.F(nF0_IDF,y0IDF_idx);

        
    end
    
    %save([dir,'\edgeSpectVals.mat'],'mEdgeSpect2','y0')
end


