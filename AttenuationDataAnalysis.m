
clear
close all 

nRepeats = 6; 
nTranslations = 10; 
omitBool = boolean(ones(1,10)); 


load('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE_StudyData\CCR5912_06V\VSX_init.mat','Trans')
bw6dB = [Trans.Bandwidth(1), Trans.Bandwidth(2)]; 

for iRepeat = 1

    frameData = load(['ABT1/Results/Frame',num2str(iRepeat),'Att.mat']);
    [~,nFm6dB]= min(abs(bw6dB(1)-frameData.BAE.IDFres.f));
    [~,nFp6dB]= min(abs(bw6dB(2)-frameData.BAE.IDFres.f));

    if 0 
        figure 
        plot(frameData.BAE.IDFres.f,frameData.BAE.IDFres.bm)
        xlabel('Frequency')
        ylabel('Normalised Attenuation (dB/cm)')
        yyaxis right 
        hold on 
        plot(frameData.BAE.IDFres.f(nFm6dB).*[1 1 ],[0 1],'k-')
        plot(frameData.BAE.IDFres.f(nFp6dB).*[1 1 ],[0 1],'k-')
        set(gca,'YColor','k')
        yticks([])
    end
    
    [slope, err, ~] = l
        

end