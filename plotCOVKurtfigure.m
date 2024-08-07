
function [] = plotCOVKurtfigure(bscEstimate,tString)
    
    %Assumes that bscEstimate has come out of COVsegmentation
    %bscEstimate = [mean(bscValues);std(bscValues);segPctLines;segPctKernels;skewness(bscValues);kurtosis(bscValues)];
     

    figure 
    sgtitle(tString)
    subplot(1,3,1)
    errorbar(bscEstimate(3,:),bscEstimate(1,:),bscEstimate(2,:),'-.')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('BSC (m^{-1}sr^{-1})')
    title('BSC')

    xlim padded
    yMax = max(bscEstimate(1,:));
    %yMin = min(bscEstimate(1,:));
    ylim padded
    tickVector = (10.^linspace(-6,log10(1.*yMax),10));
    yticks(tickVector)

    subplot(1,3,2)
    plot(bscEstimate(3,:),bscEstimate(2,:)./bscEstimate(1,:),'-x')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('C.0.V.')
    title('Coefficient of Variation')
    xlim padded
    yMax = max(bscEstimate(2,:)./bscEstimate(1,:));
    ylim padded
    %yticks(tickVector)
    
 


    subplot(1,3,3)
    plot(bscEstimate(3,:),bscEstimate(5,:),'-x')
    %set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Kurtosis')
    title('Kurtfigure')
    xlim padded
    ylim padded
    %yMax = max(COV);
    %yMin = min(COV);
    %ylim([yMin*0.5 yMax*2])
    
    
end
