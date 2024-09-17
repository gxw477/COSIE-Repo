function [iqOut,rfOut] = demodulateIQfn(PData,IQData)
    
    l=PData.Size(1);
    dz=2*PData.PDelta(3);
    nr = PData.Size(2);
    phasor=repmat([0:l-1]'*2*pi*dz,1,nr);
    phasor=reshape(phasor,l,nr);
    %IQ=interp1([0:l-1],IQData{1}(1:l,:,:),[0:0.25/dz:l-1],'cubic');

    IQa=interp1([0:l-1],abs(IQData{1}(1:l,:,:)),[0:0.25/dz:l-1],'cubic');
    IQp=interp1([0:l-1],unwrap(angle(IQData{1}(1:l,:,:))-phasor),[0:0.25/dz:l-1],'cubic');
    
    IQ=IQa.*exp(sqrt(-1)*IQp);        
    phasor=repmat(exp(sqrt(-1)*dz*2*pi*+[0:0.25/dz:l-1]'),1,nr);

    iqOut=IQa.*exp(sqrt(-1)*IQp);  


    phasor=repmat(exp(sqrt(-1)*dz*2*pi*+[0:0.25/dz:l-1]'),1,nr);
    rfOut=real(IQ.*phasor)/1896;
    

end