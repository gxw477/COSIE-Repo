
clear
close all

path
resultsDir = ['C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\SegmResults\'];


fNames = ls([resultsDir,'\SegmResults*']);
nImages = size(fNames,1)-1;

results1 = load([resultsDir,'\SegmResults1.mat']);
nEMLidx = size(results1.powerSegCOHAll{1},2)

nDepths = size(results1.allDepths,2);
allAtt_COH_all = zeros(nImages,size(results1.allAtt_COH,1),size(results1.allAtt_COH,2),size(results1.allAtt_COH,3));
allAtt_SNR_all = zeros(nImages,size(results1.allAtt_SNR,1),size(results1.allAtt_SNR,2),size(results1.allAtt_SNR,3));
allAtt_unseg = zeros(nImages,size(results1.allAtt_unseg,1),size(results1.allAtt_unseg,2),size(results1.allAtt_unseg,3));
powerSegCOHAll = zeros(nImages,nDepths,3,nEMLidx);
powerSegSNRAll = zeros(nImages,nDepths,3,nEMLidx);

for iFile = 1:nImages

    results = load([resultsDir,fNames(iFile,:)]);
    
    fields = fieldnames(results);
    
    allAtt_COH_all(iFile,:,:,:) = results.allAtt_COH;
    allAtt_SNR_all(iFile,:,:,:) = results.allAtt_SNR;

    
    
    for iDepth = 1:nDepths
        powerSegCOHAll(iFile,iDepth,1:3,:) = results.powerSegCOHAll{iDepth}([4,1,2],:);
        powerSegSNRAll(iFile,iDepth,1:3,:) = results.powerSegSNRAll{iDepth}([4,1,2],:);
    end

end


%% Plot all data 

load('C:\Users\gwest\Documents\MATLAB\EmmaLiver_NHV_NTGC\QUAD\SegmResults\CMap.mat')

mubsTH = 0.1;
mubsTH_STD = [0.04 2.2];

neg = mubsTH_STD(1);
pos = mubsTH_STD(2);

fillColor = [140, 222, 162]./256;
   
xShade = [-15 -15 110 110 ];
yShade = [ neg pos pos neg ];

for iDepth = 1:nDepths
    figure
    area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',neg,'ShowBaseLine','off'); 
    hold on
    plot([-15 110],mubsTH.*[1,1],'-.','Color',[255 171 0]./256,'LineWidth',2)
    plot([-15 110],pos.*[1 1],'-.','LineWidth',2,'Color',fillColor)
    plot([-15 110],neg.*[1 1],'-.','LineWidth',2,'Color', fillColor)
    plot([-15 -15],[neg pos],'-.','LineWidth',2,'Color',fillColor)
    plot([110 110],[neg pos],'-.','LineWidth',2,'Color',fillColor)
  
    for iFile = 1:nImages      
        errorbar(squeeze(powerSegCOHAll(iFile,iDepth,1,:)),squeeze(powerSegCOHAll(iFile,iDepth,2,:)),squeeze(powerSegCOHAll(iFile,iDepth,3,:)),'-o','Color',cmap(iFile,:),'MarkerFaceColor',cmap(iFile,:))
        hold on
        set(gca,'YScale','log')
    end
    
    h = get(gca,'Children');
    set(gca,'Children',[h(7:(7+nImages-1)); h(1:6) ])
     
    legend('1','2','3','4','5','6','7')


end


close all

%% Align all values with the same segmentation percentage

uniqueSegPct = unique(powerSegCOHAll(:,:,1,:));

powerSegCOH_pctSort = cell(1,length(uniqueSegPct));
powerSegSNR_pctSort = cell(1,length(uniqueSegPct));

for iSeg = 1:length(uniqueSegPct)
    iSeg 
    powerSegCOH_pctSort{iSeg} = [];
    powerSegSNR_pctSort{iSeg} = [];
    
    for iDepth = 1:nDepths
        
        for iFile = 1:nImages
            
            temp1 = squeeze(powerSegCOHAll(iFile,iDepth,1,:));
            temp2 = squeeze(powerSegSNRAll(iFile,iDepth,1,:));

            for iEML = 1:length(temp1)
                if temp1(iEML) == uniqueSegPct(iSeg)
                    powerSegCOH_pctSort{iSeg} = [powerSegCOH_pctSort{iSeg},squeeze(powerSegCOHAll(iFile,iDepth,1:3,iEML))] ; 
                end
                if temp2(iEML) == uniqueSegPct(iSeg)
                    powerSegSNR_pctSort{iSeg} = [powerSegSNR_pctSort{iSeg},squeeze(powerSegSNRAll(iFile,iDepth,1:3,iEML))] ; 
                end
                

            end
        end
    end
end

powerSegCOH_meanByPct = zeros(length(powerSegCOH_pctSort),3);
powerSegSNR_meanByPct = zeros(length(powerSegSNR_pctSort),3);


for iSeg = 1:length(powerSegCOH_pctSort)

    powerSegCOH_meanByPct(iSeg,1) = powerSegCOH_pctSort{iSeg}(1,1);
    powerSegCOH_meanByPct(iSeg,2) = mean(powerSegCOH_pctSort{iSeg}(2,:),'omitnan');
    powerSegCOH_meanByPct(iSeg,3) = mean(powerSegCOH_pctSort{iSeg}(3,:),'omitnan');

    if size(powerSegSNR_pctSort{iSeg},1)>0 
            powerSegSNR_meanByPct(iSeg,1) = powerSegSNR_pctSort{iSeg}(1,1);
            powerSegSNR_meanByPct(iSeg,2) = mean(powerSegSNR_pctSort{iSeg}(2,:),'omitnan');
            powerSegSNR_meanByPct(iSeg,3) = mean(powerSegSNR_pctSort{iSeg}(3,:),'omitnan');
    end
    
    
end

boolZero = powerSegSNR_meanByPct(:,1) > 0;

powerSegSNR_meanByPct = powerSegSNR_meanByPct(boolZero,:);

figure
area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',neg,'ShowBaseLine','off'); 
hold on
plot([-15 110],mubsTH.*[1,1],'-.','Color',[255 171 0]./256,'LineWidth',2)
plot([-15 110],pos.*[1 1],'-.','LineWidth',2,'Color',fillColor)
plot([-15 110],neg.*[1 1],'-.','LineWidth',2,'Color', fillColor)
plot([-15 -15],[neg pos],'-.','LineWidth',2,'Color',fillColor)
plot([110 110],[neg pos],'-.','LineWidth',2,'Color',fillColor)
errorbar(powerSegCOH_meanByPct(2:end,1),powerSegCOH_meanByPct(2:end,2),powerSegCOH_meanByPct(2:end,3),'k-o','MarkerFaceColor','k')


errorbar(powerSegSNR_meanByPct(2:end,1),powerSegSNR_meanByPct(2:end,2),powerSegSNR_meanByPct(2:end,3),'-o','Color',[0 2 255]./256,'MarkerFaceColor','blue') 
errorbar(powerSegCOH_meanByPct(1,1),powerSegCOH_meanByPct(1,2),powerSegCOH_meanByPct(1,3),'r-o','MarkerFaceColor','r')
xlabel('Final Segmentation %')
ylabel('BSC (m^{-1}sr^{-1})')

set(gca,'YScale','log')

figure 
plot(powerSegCOH_meanByPct(2:end,1),abs(powerSegCOH_meanByPct(2:end,2)- powerSegCOH_meanByPct(1,2)),'k-o','MarkerFaceColor','k')
hold on 
plot(powerSegSNR_meanByPct(2:end,1),abs(powerSegSNR_meanByPct(2:end,2)- powerSegSNR_meanByPct(1,2)),'-o','Color',[0 2 255]./256,'MarkerFaceColor','blue') 
xlabel('Final Segmentation %')
ylabel('BSC error (m^{-1}sr^{-1})')

  