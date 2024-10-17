
function [mEdgeSpect ,y0] = edgeDetectionMulti(dir,switchBool)
    % edgeDetectionMulti(dir,switchBool) 
    % dir directory 
    % switchBool = 1 for speckle directory, 0 for test directory 


    %kernel length in periods of T0 = 1/f0
    kLengthPER = 30;
    xBool = 17:111;
   
   
    IDF = load([dir,'\IDF_FiltAvg.mat']);
    
    if switchBool
        fnames = ls([dir,'\BfimgData*']);
    elseif ~switchBool
        fnames = ls([dir,'\EdgeCheck*']);
    end
    
    nImages = size(fnames,1);

    mEdgeSpect = zeros(1,nImages);
    y0 = zeros(1,nImages);
    
    for iImage = 1:nImages
        
        bfImgData = load([dir,'\',(fnames(iImage,:))]);

        fs = bfImgData.fs;
        c0 = bfImgData.Resource.Parameters.speedOfSound;
        yVals = bfImgData.yVals;
        f0 = bfImgData.Trans.frequency*1e6;
        fullIM = bfImgData.fullIM;
        
        lambda = c0/f0;
        offset = bfImgData.Receive(1).startDepth*lambda;
        wSize = IDF.Filt.wsize*lambda; 

        IDF.Filt.yVals = offset + wSize/2 + (IDF.Filt.wpos )./2*lambda/bfImgData.Receive(1).samplesPerWave;
        %wPosActual = (vsxParams.Receive(1).startDepth*sPerWave+Filt.wpos(:))/mmtosamps*1e-3; 
       
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
    
        if switchBool 
            [~,yIdx]= max(sum(fullIM,2)); 
            timeS = sum(abs(fullIM),2);
           
        else
            timeS = sum(abs(fullIM),2);
            p1 = figure ;
            plot(gca,timeS);
            hold on 
            plot(gca,sum(fullIM,2));
            xlim([100 400])
            peakMat = ginput;
            yIdx = round(mean(peakMat(:,1)));
            close(p1)
        end

        y0_Val = yVals(yIdx);

        if 1%y0_Val < min(IDF.Filt.yVals) || y0_Val > max(IDF.Filt.yVals)
            
           % mEdgeSpect(iImage) = nan;
           % y0(iImage) = nan;

        %else
            
            [~,y0IDF_idx] = min(abs(y0_Val-IDF.Filt.yVals));
            y0(iImage)= y0_Val;
            
            samplesPwavel = 4;
            kLengthSamples = kLengthPER*samplesPwavel;
            win = hanning(kLengthSamples);
        
            if yIdx-kLengthSamples/2<0 
                yIdx = kLengthSamples/2+1;
            end
        
            yBool = (yIdx - kLengthSamples/2):(yIdx+kLengthSamples/2-1);
            
            edgeRegion = fullIM(yBool,xBool).*win;
            spectVals = abs(fft(edgeRegion)).^2;
        
            df = fs/(kLengthSamples-1);
            fVals = (0:kLengthSamples-1).*df;
        
            [~ , nF] = min(abs(fVals-f0));
            
            fLine = IDF.Filt.F(nF0_IDF,:);
            idfCorrVal = interp1(IDF.Filt.yVals,fLine,yVals(yBool));

            mEdgeSpect(iImage) = mean(spectVals(nF,:))/mean(idfCorrVal,'omitnan');
            
    end
    
    %save([dir,'\edgeSpectVals.mat'],'mEdgeSpect2','y0')
end


