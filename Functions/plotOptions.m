function [] = plotOptions(data,tString,switchVble)
%plotOptions give it your data and title and it'll plot (1) kurtosis data
%or (2) coefficient of variation data
%   
    
    if switchVble == 1
        plotKurtfigure(data,tString)
    elseif switchVble ==2 
        plotCOVfigure(data,tString)
    elseif switchVble ==3
        plotCOVKurtfigure(data,tString)
    end
end