

function [] = fancySwarmPlotter(swarmAx,depth,coherence,segDepthBoolIdxs,idxsSpeckle,idxsWire,idxsCyst)
    
    c = ones(1,length(coherence));
    
    %idxsSpeckle = input('Select Speckle rejects : ');
    %idxsWire = input('Select wire rejects : ');
    %idxsCyst = input('Select cyst rejects : ');
    
    %c(idxsSpeckle) = 2; 
    %c(idxsWire) = 3; 
    %c(idxsCyst) = 4; 
    
  
    s = swarmchart(swarmAx,depth.*ones(1,length(coherence(idxsSpeckle))),coherence(idxsSpeckle),40,[43, 161, 95]./256,'filled')
    s.XJitterWidth = 2;
    hold on 
    s = swarmchart(swarmAx,depth.*ones(1,length(coherence(idxsWire))),coherence(idxsWire),40,[0,0,1],"filled")
    s.XJitterWidth = 3;
    s = swarmchart(swarmAx,depth.*ones(1,length(coherence(idxsCyst))),coherence(idxsCyst),40,[0,0,0],'filled')
    s.XJitterWidth = 4;
    s = swarmchart(swarmAx,depth.*ones(1,length(coherence(segDepthBoolIdxs))),coherence(segDepthBoolIdxs),20,[1,0,0])
    s.XJitterWidth = 1;

end
