
clear 
close all 

%Looks for which BSX estimation results are the best and tries to formally 
%state which EML point is best for all images

speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1-4Dir\';
testDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QA\';

%load verasonics param's
vsxParams = load([testDir,'\VSXoutput.mat']);
vsxParams2 = load([speckleDir,'\VSXoutput.mat']);

wOption = 2;%input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');

if wOption == 1 
    wName = 'Rect';
elseif wOption == 2
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welch';
end


% parameters
sumIdx = 33;
wName = 'Tukey';
depthSelect = 25;

resultsDirs = ls([testDir,wName,'\Z',num2str(depthSelect),'\SegResults_*']);

figure 
hold on 

ax1= subplot(2,2,1);
ax2= subplot(2,2,2);
ax3= subplot(2,2,3);
ax4= subplot(2,2,4);

hold(ax1,'on')
hold(ax2,'on')
hold(ax3,'on')
hold(ax4,'on')

for iImage = 1:3

    segDir = [testDir,wName,'\Z25\',resultsDirs(iImage,:),'\SumIdx_',num2str(sumIdx),'\'];
    segStruct = load([segDir,'SegResults.mat']);
    
    errorbar(ax1,segStruct.bscEstimate_COH_COSIE(:,3),segStruct.bscEstimate_COH_COSIE(:,1),segStruct.bscEstimate_COH_COSIE(:,2),'-.')%,'Color','k') 
    errorbar(ax2,segStruct.bscEstimate_ENV_COSIE(:,3),segStruct.bscEstimate_ENV_COSIE(:,1),segStruct.bscEstimate_ENV_COSIE(:,2),'-.')%,'Color','k')
    
    errorbar(ax3,segStruct.bscEstimate_COH_WEIGHT(:,3),segStruct.bscEstimate_COH_WEIGHT(:,1),segStruct.bscEstimate_COH_WEIGHT(:,2),'-.')%,'Color','k')
    errorbar(ax4,segStruct.bscEstimate_ENV_WEIGHT(:,3),segStruct.bscEstimate_ENV_WEIGHT(:,1),segStruct.bscEstimate_ENV_WEIGHT(:,2),'-.')%,'Color','k')
    
end


for i = 1:4
    eval(['errorbar(ax',num2str(i),',-10,segStruct.bscSpeckleBf,segStruct.bscSpeckleSTD,"r.")'])
    eval(['plot(ax',num2str(i),',[-10 100],segStruct.bscSpeckleBf.*[1 1],"r-.")'])
    eval(['xlim(ax',num2str(i),',[-15 100])'])
    
end

eval(['legend(ax',num2str(4),',"Img1","Img2","Img3","G.T")'])
    
title(ax1,'COSIE (Coh)')
title(ax2,'COSIE (SNR)')
title(ax3,'Weighting (Coh)')
title(ax4,'Weighting (SNR)')
% EML point




% vbles
% depth 
