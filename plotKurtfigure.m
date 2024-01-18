
function [] = plotKurtfigure(bscEstimate,tString)
   

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
    tickVector = round(10.^linspace(log10(10),log10(4*yMax),10));
    yticks(tickVector)

    subplot(1,3,2)
    plot(bscEstimate(3,:),bscEstimate(2,:),'-x')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Standard Deviation of BSC (a.u.)')
    title('Standard Deviation')
    xlim padded
    yMax = max(bscEstimate(2,:));
    yMin = min(bscEstimate(2,:));
    ylim([yMin*0.5 yMax*4])
    tickVector = round(10.^linspace(log10(10),log10(4*yMax),10));
    yticks(tickVector)
    
 


    subplot(1,3,3)
    plot(bscEstimate(3,:),bscEstimate(5,:),'-x')
    %set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Kurtosis')
    title('Kurtfigure')
    xlim padded
    %yMax = max(COV);
    %yMin = min(COV);
    %ylim([yMin*0.5 yMax*2])
    
    
end
