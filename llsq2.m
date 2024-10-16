function [slope,err] = llsq2(x,y)
%[slope err z] = LLSQ Computes the zero intercept of the attenuation line
    
    sqTrms = sum(x.^2);
    crossTms = sum(x.*y);

    slope = crossTms/sqTrms;
    err = sqrt(mean((y-x.*slope).^2));

    
end