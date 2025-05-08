function [] = bscSegPlotter(powerSeg,convFactor,mubsTH,mubsTH_STD,z)

    
    bscAllMean = powerSeg(1,:).*convFactor;
    bscAllStd = powerSeg(2,:).*convFactor;

    bscAllStdNeg = bscAllStd;
    
    lTHbool = bscAllMean-bscAllStd < 0;
    bscAllStdNeg(lTHbool) = bscAllMean(lTHbool);

    figure 
    errorbar(powerSeg(4,:),bscAllMean,bscAllStdNeg,bscAllStd,'k.','MarkerSize',10)
    hold on 
    plot([-15 110],mubsTH.*[1 1],'r-.','LineWidth',2)
    errorbar(-10,mubsTH,mubsTH_STD,'r.')
    set(gca,'YScale','log')
    yticks(0.125.*[1,2,4,8,16,5/0.125,10/0.125])
    xticks([0 : 25: 100])
    xlabel('Segmentation (%)')
    ylabel('BSC (m^{-1}rad^{-1})')
    set(gca,'FontSize',20)
    ylim padded
    xlim tight

    %Sorts out negative limits

    fname = ['C:\Users\gwest\Documents\COSIE paper 1\Figure_Phantom\G218L74_1ParametricImage\BSCSeg',num2str(z)];

    %exportgraphics(gcf, [fname,'.pdf'], 'ContentType', 'vector');
    %saveas(gcf,[fname,'.fig'])

end