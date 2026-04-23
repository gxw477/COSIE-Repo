
clear 
%close all 
path(path,'C:\Users\gwest\Documents\MATLAB\COSIE-Repo\Scripts\')
path(path,'C:\Users\gwest\Documents\MATLAB\COSIE-Repo\Functions\')
path(path,'C:\Users\gwest\Documents\MATLAB\AttenuationGUI\')


for volNumber = 1
    
    topDir = ['C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\'];
    
    FiltData = load([topDir,'IDFfilt.mat']);
    vsxParams = load([topDir,'\VSXoutput.mat']);
    

    analFolder = [topDir,'\AttData\'];
    %analFolder = [topDir,'\AttData\'];
    
    saveFolder =  [topDir,'\AttDataTestFolderLiverDist\'];
    %analFolder = [topDir,'/Test/'];

    if ~exist(saveFolder)
        mkdir(saveFolder)
    end
    
    frameNames = ls([analFolder,'\Frame*']);
    %frameNames = ls([analFolder,'/Att*']);
    
    nFrames = size(frameNames,1);
    
    myfittype = fittype("a0*f^b",dependent="y",independent="f",coefficients=["a0" "b"]);
    
    
    for frameIdx= 1:7  
    
        %close all
    
        load([analFolder,frameNames(frameIdx,:)])
        
        %wpos = FiltData.Filt.wpos;
        spw = 4;
        mmtosamps=2*spw/(1e-3*vsxParams.Resource.Parameters.speedOfSound/vsxParams.Trans.frequency);
    
        f= BAE.IDFres.f;


        F = BAE.IDFres.F_filt;
        
        %don't know what this function is even for 
        mup = (BAE.IDFres.F_filt(:,1));
        
        ae = median(BAE.edge_data);
        
        x = BAE.IDFres.x;

        counter = 0.1;
        
        %diffCount = (size(F,2)/1e2)-floor(size(F,2)/1e2);
        uX = unique(x);
       
        %load(['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\MRIsubcutAnalysis\MRIsegSubCut\AllDistanceMeasures.mat']);

        %slDist = DistData.allslDist(volNumber,:).*0.1; 
        
        figure 
        iq = load([topDir,'/BFimgData',num2str(frameIdx),'.mat'],'iq','yVals');
        b80 = bmode(iq.iq,80);
        imagesc(1:128,iq.yVals,b80')
        colormap gray

        pause(0.1)
        
        slDist = input('Depth (cm) : ');

        nDepth = length(uX(:));
        [~, depthMinU ] = min(abs(slDist-uX));       
        
        nPerDepth = floor(length(x)/length(uX));
        
        ixStart = (depthMinU-1)*nPerDepth+1;


        xEst = zeros(1,nDepth-1);
        m0dB = zeros(1,nDepth-1);
        m0dBErr = zeros(1,nDepth-1);
        aX = zeros(nDepth-1,length(f));
        bX = zeros(nDepth-1,length(f));
        eX= zeros(nDepth-1,length(f),2);
        re = zeros(nDepth-1,length(f));
    
        [~,fL] = min(abs(4- f));
        [~,fU] = min(abs(6.5 - f));
        
        a0 = zeros(1,nDepth-1);
        alpha = a0;
        corr = a0;
        err =zeros(nDepth-1,2);
        res = zeros(1,nDepth-1);
        
        nValues = nPerDepth*nDepth;
        
        counter = 0.1; 
        
        for iX = 2:nDepth
        
            
            iXEnd = (iX -1) * nPerDepth ;

            a = zeros(1,length(f));
            b = a;
            c = a; 
            e = zeros(length(f),2);
            re = a;
            
            %Attenuation calculation
            
            for iF=1:length(f)         
                
              
                [a(iF),b(iF),c(iF),e(iF,:),re(iF)]=lsqn(2*x(ixStart:iXEnd),10*log10(F(iF,ixStart:iXEnd)),[BAE.ABTpos 0 0]);
    
                if iX*nPerDepth == size(F,2)
                    figure
                    plot(2*x(1:iXEnd),10*log10(F(iF,1:iXEnd)),'x')
                    hold on
                    plot(2*x(1:iXEnd),a(iF) + 2*x(1:iXEnd).*b(iF),'r-')
                    close all
                end
                    
            end
    
            xEst(iX-1) = x(iXEnd);
            aX(iX-1,:) = a;
            bX(iX-1,:) = -b;
            eX(iX-1,:,:) = e;
            re(iX-1,:) = re;
    
            
            [a0(iX-1),alpha(iX-1),corr(iX-1)] = lsqn(f(fL:fU),bX(iX-1,fL:fU),[f(fL) 0 0]); 
            
            temperr = (bX(iX-1,fL:fU)-(a0(iX-1)+alpha(iX-1))).^2;
            err(iX-1) = sqrt(sum(temperr));


            m0dB(iX-1) = sum(f(fL:fU).*bX(iX-1,fL:fU))./sum(f(fL:fU).^2);
            
            y0dB = f(fL:fU).*m0dB(iX-1);

            m0dBErr(iX-1) = sqrt(sum((bX(iX-1,fL:fU) -y0dB).^2));
            
            
            
            if iX > nDepth*counter
                x(iXEnd-1:iXEnd+1);
                [volNumber, frameIdx,round(iX/nDepth*100)]
                counter = counter +  0.1;
                
                figure 
                plot(f(fL:fU),bX(iX-1,fL:fU),'k.')
                hold on
                plot(f(fL:fU),f(fL:fU).*m0dB(iX-1),'r-.')
                plot(f(fL:fU),f(fL:fU).*alpha(iX-1)+a0(iX-1),'r-*')
                
                close all
               

            end

        end
    
    

        figure 
        plot(sort(xEst),m0dB,'k.')
        hold on 
        plot(sort(xEst),alpha,'k*')
        xlabel('End Depth of attenuation estimate (cm)')
        ylabel('Attenuation slope f.t. 0 (dB/cm/MHz)')

        
        figure 
        plot(sort(xEst),m0dB.*vsxParams.Trans.frequency,'k.')
        hold on 
        plot(sort(xEst),alpha.*vsxParams.Trans.frequency + a0,'k*')
        xlabel('End Depth of attenuation estimate (cm)')
        ylabel('Attenuation p.u. distance (dB/cm)')




        close all

        save([saveFolder,'/attFit',num2str(frameIdx),'.mat'],'xEst','alpha','a0','f','err','aX','bX','eX','re','m0dB','m0dBErr')
    end
end
