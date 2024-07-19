
clear
close all

topDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\QA\';%[uigetdir('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74\','Select Analysis directory'),'\'];

%speckleDir = ['C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\PhantomExperimentsL74\Speckle\'];
speckleDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\ElastPhtL74_1607\Img1Dir\';

%topDir = speckleDir;

sumIdx = 32;
depthIdx = 25;% input('Depth (mm) : '); 

iImage = input('Image # :' );
wOption = input('Window Type \n 1 for rectangular \n 2 for tukey \n 3 for Welch : \n ');

if wOption == 1 
    wName = 'Rect';
elseif wOption == 2
    wName = 'Tukey';
elseif wOption == 3
    wName = 'Welch';
end


nImageSpeckle = 5; 

allSpeckleCOV = zeros(1,nImageSpeckle);
allSpeckleKurt = zeros(1,nImageSpeckle);

for iImageSpeckle = 1:nImageSpeckle

    load([speckleDir,wName,'\Z',num2str(depthIdx),'\SegResults_',num2str(iImageSpeckle),'_adaptive\SumIdx_',num2str(sumIdx),'\StatAttack.mat'])

    allSpeckleCOV(iImageSpeckle) = bscEstimate_COH(2,1)./bscEstimate_COH(1,1);
    allSpeckleKurt(iImageSpeckle)= bscEstimate_COH(5,1);

end


speckleCOV_stats = [mean(allSpeckleCOV), std(allSpeckleCOV)];
speckleKurt_stats = [mean(allSpeckleKurt), std(allSpeckleKurt)];

load([topDir,wName,'\Z',num2str(depthIdx),'\SegResults_',num2str(iImage),'_adaptive\SumIdx_',num2str(sumIdx),'\StatAttack'])


figure 
subplot(1,2,1)
plot(bscEstimate_COH(3,:),bscEstimate_COH(2,:)./bscEstimate_COH(1,:),'-o','Color','k','MarkerFaceColor','k')
hold on 
plot(bscEstimate_weighted_COH(3,:),bscEstimate_weighted_COH(2,:)./bscEstimate_weighted_COH(1,:),'-sq','Color','k','MarkerFaceColor','k')
plot(bscEstimate_ENV(3,:),bscEstimate_ENV(2,:)./bscEstimate_ENV(1,:),'-o','Color','red','MarkerFaceColor','red')
plot(bscEstimate_weighted_ENV(3,:),bscEstimate_weighted_ENV(2,:)./bscEstimate_weighted_ENV(1,:),'-sq','Color','red','MarkerFaceColor','red')

plot([0 100],speckleCOV_stats(1).*[1 1],'-.','Color','magenta')
errorbar(100,speckleCOV_stats(1),speckleCOV_stats(2),'o','Color','magenta','MarkerFaceColor','auto')
xlim padded
xlabel('Seg %')
ylabel('C.O.V of BSC')
legend({'COSIE-COH','Weighted-COH','COSIE-ENV','Weighted-ENV'})

ax = gca; 
ax.Legend.Position =  [0.5783 0.7448 0.1330 0.1655];


subplot(1,2,2)
plot(bscEstimate_COH(3,:),bscEstimate_COH(5,:),'-o','Color','k','MarkerFaceColor','k')
hold on 
plot(bscEstimate_weighted_COH(3,:),bscEstimate_weighted_COH(5,:),'-sq','Color','k','MarkerFaceColor','k')
plot(bscEstimate_ENV(3,:),bscEstimate_ENV(5,:),'-o','Color','red','MarkerFaceColor','red')
plot(bscEstimate_weighted_ENV(3,:),bscEstimate_weighted_ENV(5,:),'-sq','Color','red','MarkerFaceColor','red')

plot([0 100],speckleKurt_stats(1).*[1 1],'-.','Color','magenta')
errorbar(100,speckleKurt_stats(1),speckleKurt_stats(2),'o','Color','red','MarkerFaceColor','auto','MarkerSize',6)
xlim padded
xlabel('Seg %')
ylabel('Kurtosis of BSC')
legend({'COSIE-COH','Weighted-COH','COSIE-ENV','Weighted-ENV'})

ax = gca; 
ax.Legend.Position =  [0.5783 0.7448 0.1330 0.1655];


%saveas(gcf,['QAphtAngleCKfigs/img',num2str(iImage2)],'png')

