function [slope,err,z] = llsq(BAE,flag)
%[slope err z] = LLSQ Computes the zero intercept of the attenuation line
%  slope = sum(x.*y)/sum(x.^2)
%  err = y - x*slope
%  z = depth of ABT window
%  flag == 1 : ABT
%  flag == 2 : IDF  
    
    if flag == 1 
        x = BAE.ABTres.fs;
        y = BAE.ABTres.bs;
    elseif flag ==2
        x = BAE.IDFres.fs;
        y = BAE.IDFres.bs;
    end
    
    sqTrms = sum(x.^2);
    crossTms = sum(x.*y);

    slope = crossTms/sqTrms;
    err = sqrt(mean((y-x.*slope).^2));
    z = BAE.ABTpos;
    
end