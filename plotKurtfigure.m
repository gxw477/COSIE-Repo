
function [] = plotKurtfigure(bscEstimate,tString)
   

    figure 
    sgtitle(tString)
    subplot(1,3,1)
    errorbar(bscEstimate(3,:),bscEstimate(1,:),bscEstimate(2,:),'-.')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Backscattered power (a.u.)')
    title('Power')
    xlim([-5 max(bscEstimate(3,:))+10])
    yMax = max(bscEstimate(1,:));
    yMin = min(bscEstimate(1,:));
    ylim([yMin*0.5 yMax*4])
    tickVector = round(10.^linspace(log10(10),log10(4*yMax),10));
    yticks(tickVector)

    subplot(1,3,2)
    plot(bscEstimate(3,:),bscEstimate(2,:),'-x')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Standard Deviation of BSC (a.u.)')
    title('Standard Deviation')
    xlim([-5 max(bscEstimate(3,:))+10])
    yMax = max(bscEstimate(2,:));
    yMin = min(bscEstimate(2,:));
    ylim([yMin*0.5 yMax*4])
    tickVector = round(10.^linspace(log10(10),log10(4*yMax),10));
    yticks(tickVector)
    
 


    subplot(1,3,3)
    COV = bscEstimate(2,:)./bscEstimate(1,:);
    plot(bscEstimate(3,:),COV,'-x')
    set(gca,'YScale','log')
    xlabel('Segmentation %')
    ylabel('Kurtosis')
    title('Kurtfigure')
    xlim([-5 max(bscEstimate(3,:))+10])
    yMax = max(COV);
    yMin = min(COV);
    ylim([yMin*0.5 yMax*2])
    


    tickVector = linspace(0,3,12);
    yticks(tickVector)

end
