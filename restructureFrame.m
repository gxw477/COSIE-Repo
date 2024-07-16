

function [frameOUT ] = restructureFrame(frameIN,connections)

    frameOUT = zeros(size(frameIN));
    
    for iChannel = 1:length(connections)
        
        channelOut = connections(iChannel);
        frameOUT(:,iChannel) = frameIN(:,channelOut);

    end


end