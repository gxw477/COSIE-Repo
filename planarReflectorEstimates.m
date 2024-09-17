
function R0 = planarReflectorEstimates(dir,prefDir,switchBool)

    
    kLengthPER = 30;
    xBool = 17:111;
    
    fnames = ls([dir,'\EdgeCheck*']);
    
    nImages = size(fnames,1);

    mEdgeSpect = zeros(1,nImages);
    y0 = zeros(1,nImages);
    
    dVal = [10 15]; 


    for iImage = 1:length(dVal)
        
        bfImgData = load([dir,'\',(fnames(iImage,:))]);

        fs = bfImgData.fs;
        c0 = bfImgData.Resource.Parameters.speedOfSound;
        yVals = bfImgData.yVals;
        f0 = bfImgData.Trans.frequency*1e6;
        fullIM = bfImgData.fullIM;    
        
        lambda = c0/f0;
        offset = bfImgData.Receive(1).startDepth*lambda;
        
        if size(fullIM,2)>size(fullIM,1)
            fullIM = fullIM';
        end
        

        if switchBool 
            [~,yIdx]= max(sum(fullIM,2)); 
            timeS = sum(abs(fullIM),2);
           
        else
            timeS = sum(abs(fullIM),2);
            p1 = figure ;
            plot(gca,timeS);
            hold on 
            plot(gca,sum(fullIM,2));
            xlim([100*iImage 500*iImage])
            peakMat = ginput;
            yIdx = round(mean(peakMat(:,1)));
            close(p1)
        end

        y0_Val = yVals(yIdx);
            
        prefImData = load([prefDir,'BFimgData',num2str(dVal(iImage)),'.mat']);
        prefImg = prefImData.fullIM';

        
        samplesPwavel = 4;
        kLengthSamples = kLengthPER*samplesPwavel;
        win = hanning(kLengthSamples);


        df = fs/(kLengthSamples-1);
        fVals = (0:kLengthSamples-1).*df;
        [~ , nF] = min(abs(fVals-f0));
      
        
        yBool = (yIdx - kLengthSamples/2):(yIdx+kLengthSamples/2-1);
        
        edgeRegionT = fullIM(yBool,xBool);%.*win;
        spectValsT = abs(fft(edgeRegionT)).^2;

        edgeRegionP = prefImg(yBool,xBool);%.*win;
        spectValsP = abs(fft(edgeRegionP)).^2;
        
        perspWint = 0.1307;

        R0(iImage) = mean(spectValsT(nF,:))./mean(spectValsP(nF,:))* perspWint; 

        
    end
        
   
end