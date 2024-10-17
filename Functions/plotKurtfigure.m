
function [] = plotKurtfigure(bscEstimate,tString)
    
    %Assumes that bscEstimate has come out of COVsegmentation
    %bscEstimate = [mean(bscValues);std(bscValues);segPctLines;segPctKernels;skewness(bscValues);kurtosis(bscValues)];
     

    figure 
    sgtitle(tString)
    subplot(1,3,1)
    errorbar(bscEstimate(3,:),bscEstimate(1,:),bscEstimate(2,:),'-.')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Backscattered power (a.u.)')
    title('Power')
    xlim padded
    yMax = max(bscEstimate(1,:));
    %yMin = min(bscEstimate(1,:));
    ylim padded
    tickVector = (10.^linspace(-6,log10(1.*yMax),10));
    yticks(tickVector)

    subplot(1,3,2)
    plot(bscEstimate(3,:),bscEstimate(2,:),'-x')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Standard Deviation of BSC (a.u.)')
    title('Standard Deviation')
    xlim padded
    yMax = max(bscEstimate(2,:));
    ylim padded
    tickVector = (10.^linspace(-6,log10(1.*yMax),10));
    yticks(tickVector)
    
 


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
