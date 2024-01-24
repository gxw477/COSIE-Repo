

%loads beamformed image data and select ray lines and depths of features 
clear 
close all 

toptopDir = ['PhantomExperimentsL74_QuadInterp/']

contBool = 1; 

while contBool 

    [imgFName,imgDirName] = uigetfile(toptopDir);
    fullImgName = [imgDirName,imgFName];

    path(path,imgDirName)
    imgData = load(fullImgName);
    B80 = bmode(imgData.iq',80);

    %lWidth = (imgData.Trans.ElementPos(2,1)- imgData.Trans.ElementPos(1,1))*imgData.lambda;

    imagesc(B80);
    xlabel('Lateral Position (ray)')
    axis tight
    ylabel('Axial Position (sample)')
    pbaspect([0.5,1,0.5])
    
    pause(1e-3)
    
    contBool2 = 1;
    
    while contBool2 
    
        [xFeature, yFeature] = ginput;
        
        if length(xFeature) == 2*round(length(xFeature)/2)
            contBool2 =0;
        end
    end
    

    xFeature2 = [xFeature(1:2:end -1 ), xFeature(2:2:end)];
    yFeature2 = 0.5.*(yFeature(1:2:end-1) + yFeature(2:2:end));
    
    save([imgDirName,'\Features',imgFName(length(imgFName)-4)],'xFeature2','yFeature2')

    %contBool = input('Continue ? :');

end