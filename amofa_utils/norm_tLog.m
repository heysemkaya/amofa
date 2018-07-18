function tLogN=norm_tLog(tLog,numSamples,numMeans)

    if size(tLog,2)==1
        tLogN = ones(numSamples,1);
    else
        tLogN = tLog ./ repmat((sum(tLog'))',1,size(tLog,2));
        if sum(sum(isnan(tLogN)))>0
            %somewhere, there was a division by zero...
            dummyProb = 1/numMeans; %uniform prob...
            for i=1:size(tLogN,1)
                if isnan(tLogN(i,1))
                    % once one is NaN, the whole line is NaN
                    tLogN(i,:) = ones(1,numMeans)*dummyProb;
                end
            end
        end
        if min(sum(tLogN'))<0.9999
            %there is a too small log somewhere...
            smallLogIndices = find((sum(tLogN,2)<0.1));
            for i=smallLogIndices
                [dummy maxEl] = max(tLog(i,:));
                tLogN(i,maxEl)=1;
                tLogN(i,:) = tLogN(i,:) / sum(tLogN(i,:));
            end
        end
    end %normalize tLog
            
end

