function [] = bscSegPlotter(powerSeg,convFactor,mubsTH,mubsTH_STD,z)

    
    [~,uIDs,~] = unique(powerSeg(4,:));

    powerSeg = powerSeg(:,uIDs);
    
    fillColor = [140, 222, 162]./256;

    bscAllMean = powerSeg(1,:).*convFactor;
    bscAllStd = powerSeg(2,:).*convFactor;

    bscAllStdNeg = bscAllStd;
    
    lTHbool = bscAllMean-bscAllStd < 0;
    bscAllStdNeg(lTHbool) = bscAllMean(lTHbool);


    if length(mubsTH_STD)==1
        mubsTHneg = - mubsTH_STD;
        mubsTHpos = mubsTH_STD;
    elseif length(mubsTH_STD)==2
        mubsTHneg = mubsTH_STD(1);
        mubsTHpos = mubsTH_STD(2);
    end    

    figure 
    errorbar(powerSeg(4,1),bscAllMean(1),bscAllStdNeg(1),bscAllStd(1),'r.','MarkerSize',10) 
    hold on
    errorbar(powerSeg(4,2:end),bscAllMean(2:end),bscAllStdNeg(2:end),bscAllStd(2:end),'k.','MarkerSize',10) 
    plot([-15 110],mubsTH.*[1 1]+ mubsTHpos,'-.','LineWidth',2,'Color',fillColor)
    plot([-15 110],mubsTH.*[1 1]- mubsTHneg,'-.','LineWidth',2,'Color', fillColor)
    plot([-15 -15],mubsTH + [-mubsTHneg mubsTHpos],'-.','LineWidth',2,'Color',fillColor)
    plot([110 110],mubsTH + [-mubsTHneg mubsTHpos],'-.','LineWidth',2,'Color',fillColor)
    plot([-15 110],mubsTH.*[1,1],'-.','Color','b','LineWidth',2)
    

    xShade = [-15 -15 110 110 ];
    yShade = mubsTH + [- mubsTHneg mubsTHpos mubsTHpos -mubsTHneg ];
    
   
    
    area(xShade, yShade,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',fillColor,'BaseValue',mubsTH-mubsTHneg,'ShowBaseLine','off');

    

    
    %errorbar(-10,mubsTH,mubsTH_STD,'r.')
    set(gca,'YScale','log')
    
    %yticks(0.125.*[1,2,4,8,16,5/0.125,10/0.125])
    %yticks(0.05.*2.^[-1:6])

    xticks([0 : 25: 100])
    xlabel('Segmentation (%)')
    ylabel('BSC (m^{-1}sr^{-1})')
    set(gca,'FontSize',20)
    ylim padded
    xlim([-14 104])


    %Sorts out negative limits

    fname = ['C:\Users\gwest\Documents\COSIE paper 1\Figure_Phantom\G218L74_1ParametricImage\BSCSeg',num2str(z)];

    %exportgraphics(gcf, [fname,'.pdf'], 'ContentType', 'vector');
    %saveas(gcf,[fname,'.fig'])

end

