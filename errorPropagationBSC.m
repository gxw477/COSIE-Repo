function [errs ] = errorPropagationBSC(Si,eSi,Sp,eSp,dzT,ealphaT,dzR,ealphaR,bscR,ebscR,mubsT)

    %1 error value per quantity (all variables plus final answer)
    errs = zeros(1,6);

    perrSi = (mubsT/Si)*eSi;
    perrSp = (mubsT/Sp)*eSp;
    perralphaT = (mubsT*dzT)*ealphaT*log(10);
    perralphaR = (mubsT*dzR)*ealphaR*log(10);
    perrbscR = (mubsT/bscR)*ebscR;
    
    errs = [perrSi,perrSp,perralphaT,perralphaR,perrbscR];

    errs(6) = sqrt(sum(errs(1:5).^2));
    

end