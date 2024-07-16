

%loads beamformed image data and select ray lines and depths of features 
clear 
close all 

toptopDir = ['PhantomExperimentsL74_QuadInterp/']


contBool = 1; 

while contBool 

    [imgFName,imgDirName] = uigetfile(toptopDir);
    fullImgName = [imgDirName,imgFName];

    vsxParams = load([imgDirName,'/VSXoutput.mat']);

    path(path,imgDirName)
    imgData = load(fullImgName);
    B80 = bmode(imgData.iq',80);

    %lWidth = (imgData.Trans.ElementPos(2,1)- imgData.Trans.ElementPos(1,1))*imgData.lambda;

    imagesc(B80);
    xlabel('Lateral Position (ray)')
    axis tight
    ylabel('Axial Position (sample)')
    pbaspect([0.5,1,0.5])
    colormap gray

    pause(1e-3)

    contBool2 = 1;
    mask = zeros(size(B80));

    while contBool2 
        
        m  = imfreehand('Closed',1);
        maskCurrent = m.createMask;
        
        mask = mask | maskCurrent;

        contBool2 = input('Continue? : 1 for yes, 0 for no ');

    end
    
    close all

    lWidth = (vsxParams.Trans.ElementPos(2,1)-vsxParams.Trans.ElementPos(1,1))*vsxParams.lambda;
    xVals = (1:vsxParams.P.numRays).*lWidth; 
    yVals = vsxParams.lambda.*linspace(vsxParams.Receive(1).startDepth,vsxParams.Receive(1).endDepth,size(B80,1));%
   
    
    [xQ , yQ] = meshgrid(xVals,yVals);

    
    
    save([imgDirName,'\Features',imgFName(length(imgFName)-4)],'mask')

    contBool = input('Continue ? :');

end