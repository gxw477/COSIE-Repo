function [] = EMLvideoMaker(saveDirJPGS,nEMLpoints,cohTest,speckleCOSIE,lWidth,yVals,bfImgData,depthIdx,axIdxs,powf0,xBool,sumIdx)
    
    if ~exist(saveDirJPGS,'dir')
        mkdir(saveDirJPGS)
    end
    
    v = VideoWriter([saveDirJPGS,'myFile.mp4']);
    v.FrameRate = 1.5;
    open(v)
    
    for idxEML = 1:nEMLpoints

        segBool3 = cohTest > speckleCOSIE.redEML(1,idxEML) & cohTest < speckleCOSIE.redEML(2,idxEML);
        bmCohCOSIEparImage_EML_video_FRAME((1:128).*lWidth,yVals,bfImgData,depthIdx,axIdxs,powf0,segBool3,xBool,speckleCOSIE,cohTest,sumIdx,idxEML)
        saveas(gcf,[saveDirJPGS,'idx_',num2str(idxEML),'.jpg']);
        A = imread([saveDirJPGS,'idx_',num2str(idxEML),'.jpg']);
        writeVideo(v,A);
        close all
    
    end
    
    close(v)

end