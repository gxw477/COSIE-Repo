
clear
close all 

nRepeats = 5; 
nTranslations = 6; 
omitBool = boolean(ones(1,10)); 
topDir = 'C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE_StudyData\C16D_ABT\';

load('C:\Users\gwest\Documents\Vantage-4.9.2-2308102000\COSIE_StudyData\CCR5912_06V\VSX_init.mat','Trans')
bw6dB = [Trans.Bandwidth(1), Trans.Bandwidth(2)]; 


attSpeckle  = [0.524 , 0.09];

attData = zeros(nRepeats,nTranslations);
errData = zeros(nRepeats,nTranslations);

for iRepeat = 1:nRepeats
        
    for iTranslations = 1:6
        %frameData = load([topDir,'ABT',num2str(iRepeat),'/Results/Frame',num2str(iTranslations),'Att.mat']);
        frameData = load([topDir,'ABT',num2str(iRepeat),'/Results/AttFrame',num2str(iTranslations),'.mat']);
        
        [~,nFm6dB]= min(abs(bw6dB(1)-frameData.BAE.IDFres.f));
        [~,nFp6dB]= min(abs(bw6dB(2)-frameData.BAE.IDFres.f));
        freqVector = nFm6dB:nFp6dB;

        f = frameData.BAE.IDFres.f(freqVector);
        b = frameData.BAE.IDFres.bm(freqVector);
        e = frameData.BAE.IDFres.c(freqVector);
        [slope, err ] = llsq2(f, b);
        
        attData(iRepeat,iTranslations) = slope;
        errData(iRepeat,iTranslations) = err;


        if 0 
            figure 
            errorbar(f,b,e)
            xlabel('Frequency')
            ylabel('Attenuation (dB/cm)')
            hold on 
            plot([0 , frameData.BAE.IDFres.f(freqVector(1))] , [0 , frameData.BAE.IDFres.f(freqVector(1))].*slope ,'r-.')
            plot([frameData.BAE.IDFres.f(freqVector(1)),frameData.BAE.IDFres.f(freqVector(end))] ,[frameData.BAE.IDFres.f(freqVector(1)),frameData.BAE.IDFres.f(freqVector(end))].*slope ,'r-')
    
            yyaxis right 
            hold on 
            plot(frameData.BAE.IDFres.f(nFm6dB).*[1 1 ],[0 1],'k-')
            plot(frameData.BAE.IDFres.f(nFp6dB).*[1 1 ],[0 1],'k-')
            set(gca,'YColor','k')
            yticks([])
            set(gca,'FontSize',15)
        end

    end    
    
end



%% Depth Dependence 

depth = 0.1.*(10:5:35);

% all data 
figure
hold on 
for iRepeat = 1:nRepeats
    errorbar(depth,attData(iRepeat,:),errData(iRepeat,:),'o','MarkerSize',3)
end
ylabel('Attenuation coeff (dB/cm/MHz)')
xlabel('Edge Depth (cm)')
set(gca,'FontSize',15)
xlim padded


%mean 
meanAttDepth = mean(attData,1);
meanErrDepth = mean(errData,1);

[slope,unc] = llsq2(depth,meanAttDepth);
P = polyfit(depth,meanAttDepth,1);


figure 
errorbar(depth,meanAttDepth,meanErrDepth,'ko','MarkerSize',3)
hold on 
errorbar(0,attSpeckle(1),attSpeckle(2),'ro','MarkerSize',3)
plot([0 max(depth)+0.5],attSpeckle(1).*[1 1],'r-.')
plot(depth,P(2)+depth.*P(1),'r-')
legend({'Mean \alpha',[num2str(round(P(2),3)),' ',num2str(round(P(1),3)),'*z']})
set(gca,'FontSize',15)
xlim padded 
ylabel('\alpha (dB/cm/MHz)')
xlabel('Edge Depth (cm)')

%% Repeat inDependence 


rVals = 1:nRepeats;
meanAttRep = mean(attData,2);
meanErrRep = mean(errData,2);

figure 
errorbar(rVals,meanAttRep,meanErrRep,'ko','MarkerSize',3)
hold on 
errorbar(0,attSpeckle(1),attSpeckle(2),'ro','MarkerSize',3)
plot([-1 max(rVals)+1],attSpeckle(1).*[1 1],'r-.')
xlim tight
ylabel('\alpha (dB/cm/MHz)')
set(gca,'FontSize',15)
xticks(0:1:6)
xticklabels({'G.T','1','2','3','4','5','6'})

