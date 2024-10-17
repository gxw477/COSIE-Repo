function [R] = CoherenceAnalysisFN(scatData)

%data -> (nSamples, nRec)

    nSamples = size(scatData,1);
    nRec = size(scatData,2);

    if 0% nRec > nSamples 
        scatData = scatData';

        nSamples = size(scatData,1);
        nRec = size(scatData,2);

    end
    
    mProduct = zeros(nRec,nSamples,nRec);


     for iRec1 = 1:nRec 

                mProduct(:,:,iRec1) =  (scatData.*scatData(:,iRec1))';

     end

    sumMProduct = sum(mProduct,2);

    %reshape into an nRec x nRec grid
    sigma2Terms = reshape(sumMProduct,[nRec,nRec]);
    
    %Squared terms are the diagonals 
    sqTerms = diag(sigma2Terms);

    mVals = 1:(nRec);

    R = zeros(nRec,1);

    for iLag = 1:length(mVals)-1


        term4 = zeros(1,nRec - mVals(iLag));

        for iNode = 1:nRec-mVals(iLag)

            term1 = sigma2Terms(iNode,iNode+mVals(iLag));
            term2 = sqTerms(iNode);
            term3 = sqTerms(iNode+mVals(iLag));

            if term2 == 0 || term3 == 0
                term4(iNode) = 0;
            else
                term4(iNode) = term1/(term2*term3)^0.5;
            end
        end

        R(iLag) = (1/(nRec-mVals(iLag)))*sum(term4);
        term4 = term4.*0;

    end

end