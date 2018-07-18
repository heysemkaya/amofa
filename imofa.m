function [mixture,trainingHistory] =imofa(data, cv)
%[mixture,trainingHistory]=imofa(training, validation);
% mixture has fields Lambda, Psi, Mu, Pi, numFactors
% trainingHistory has fields logsTra, logsCv, progress, totalEM
% numFactors initially indicates minimum number of factors per component
% Albert Ali Salah
% 1-factor analytic solution based on Z. Ghahramani's ffa code

cp = 3; %Mardia multivariate kurtosis based component selection
regTerm = 1; %regularization term for determinant problems
numFactors = 1; %initial number of factors per component
maxIter = 100; %max. EM iterations per run
maxTrial = 3; %try cont_mfa 3 times before giving up

numDims = size(data,2);
numSamples = size(data,1);
numMeans = 1;
smallestComponent = numDims;
lambdaSmall = 0.001;
dd = numDims * (numDims+2); %for Mardia kurtosis metric
totalEM = 0; %the total number of EM iterations

%fit a 1-factor, 1-component model
ns = '1';
eval(['[Lambda',ns,',Psi',ns,',logsTra] = ffa3(data,numFactors,maxIter);']);
eval(['Mu',ns,'=mean(data);']);
Pi=1;
totalEM = totalEM + size(logsTra,2); %add to total number of iterations

Mu1 = Mu1'; %means are stored as column vectors

oldCvLog = -Inf;
CvLog = loglike(cv,Lambda1,Psi1,Mu1); %the last loglike for validation set
logsCv = CvLog; %keeps track of validation loglike during the whole algorithm
tLog = exp(loglike2(data,Lambda1,Psi1,Mu1)); %the total weighted likelihoods for each sample
%normalize tLog
tLogN = ones(numSamples,1);%normalized posteriors

%create the initial batch variables
Lambda = Lambda1;
Psi = Psi1;
Mu = Mu1;

progress(1,1) = size(logsTra,2); %1st column keeps track of iterations
progress(1,2) = 0; %which component
progress(1,3) = 0; %what is done
progressCount = 2;

while CvLog > oldCvLog
    oldCvLog = CvLog;
    %look at the component with lowest average kurtosis for splitting
    % and largest Psi values for factor addition
    factorCount = 0;
    if numMeans == 1
        componentToSplit = 1;
        componentToAddFactor = 1;
        if numFactors(1)>numDims*.8
            componentToAddFactor = -1;
        end
    else
        if cp==2
            %soft geometric mean likelihood based split
            splitMetric = sgml2(tLog,tLogN,Lambda,Psi,numFactors);
        end

        factorCount2 = 0;
        for k=1:numMeans
            ks = num2str(k);
            %find the data set that is associated with this component
            dat = zeros(1,numDims);
            datCount = 1;
            for i=1:numSamples
                [tx ty] = max(tLog(i,:));
                if ty==k
                    dat(datCount,:) = data(i,:);
                    datCount = datCount + 1;
                end
            end
            kD=factorCount2+1:factorCount2+numFactors(k);
            factorCount2 = factorCount2+numFactors(k);
            tLambda = Lambda(:,kD);
            tempCov = tLambda*tLambda' +diag(Psi(:,k));

            if cp==1
                %univariate kurtosis metric
                splitMetric(k) = sum(abs(kurto(dat) - 3));
            elseif cp==3
                %Mardia multivariate kurtosis metric
                invCov = inv(tempCov);
                tMu = Mu(:,k)';
                b_sum = 0;
                for i=1:numSamples
                    x = data(i,:) - tMu;
                    b_sum = b_sum + tLogN(i,k)*(x*invCov*x')^2;
                end
                b_sum = b_sum / sum(tLogN(:,k));
                splitMetric(k) = abs((b_sum - dd) / sqrt(8*dd/sum(tLogN(:,k))));
            end

            % factor metric: the sum of squares of the covariance differences
            factorMetric(k) = sum(sum(((tempCov - cov(dat)).*(1-eye(size(tempCov,1)))).^2));

            %if the size of the component is too small, do not split it
            if Pi(k)*numSamples <smallestComponent
                splitMetric(k) = -Inf;
            end
            %if the number of factors are too much, do not consider adding new factors
            if numFactors(k)>numDims*.8
                factorMetric(k) = -Inf;
            end
        end
        [dummy componentToSplit] = max(splitMetric);
        %greatest total deviation from 3 in kurtosis or greatest deviation from zero
        [dummy componentToAddFactor] = max(factorMetric); %greatest covariance difference
        if max(factorMetric)==-Inf
            %we don't want to add factors
            componentToAddFactor = -1;
        end
        if max(splitMetric)==-Inf
            %we don't want to add components
            componentToSplit = -1;
        end
    end

    %initialize newLogs, to select which action pays off better
    newLogs = zeros(1,2);

    %split the componentToSplit
    k = componentToSplit;
    if k>-1
        ks = num2str(k);
        %first find the split
        %eval(['[X1 X2 split_g init_log] = softsplit(data,Lambda',ks,',Psi',ks,',tLogN(:,k));']);
        [ll indic]=max(tLogN,[],2);
        cluster_k=(indic==k);
        localdata=data(cluster_k,:);
        eval(['[X1 X2 split_g init_log] = softsplit(localdata,Lambda',ks,',Psi',ks,',tLogN(cluster_k,k),Mu',ks,');']);
        %create a dummy set of batch variables
        dumLambda = [Lambda Lambda(:,factorCount+1:factorCount+numFactors(k))];
        dumNumFactorsSplit = [numFactors numFactors(k)];
        dumPsi = [Psi Psi(:,k)];
        dumMu = [Mu X2];
        dumMu(:,k) = X1;
        dumPi = Pi;
        dumPi(k,1) = dumPi(k,1) *.5;  %the prior of the component is halved
        dumPi(numMeans+1,1) = dumPi(k,1);
        %alternatively, we can use split_g
        %dumPi(numMeans+1,1) = dumPi(k,1) * split_g(2);
        %dumPi(k,1) = dumPi(k,1) * split_g(1);

        trials = 0;
        found = 0;
        while (trials<maxTrial)
            try
                eval(['[dumLambdaSplit,dumPsiSplit,dumMuSplit,dumPiSplit,dumLogsSplit] = cont_mfa2(data,dumMu,dumLambda,dumPsi,dumPi,dumNumFactorsSplit,init_log,maxIter);']);
                totalEM = totalEM + size(dumLogsSplit,2); %add to total number of iterations

                %calculate the likelihood of this new model on the cv set
                newLogs(1) = loglike3(cv,dumLambdaSplit,dumPsiSplit,dumMuSplit,dumPiSplit,dumNumFactorsSplit);
                %set trials to maxTrial to exit while loop
                trials = maxTrial;
                found = 1; %we had no problems
            catch
                %increment the trial counter
                disp('Warning: cont_mfa2 could not split');
                trials = trials + 1;
            end
        end
        %check whether the trials succeeded
        if ~found
            newLogs(1) = -Inf;
        end
    else
        newLogs(1) = -Inf;
    end

    %now try adding a factor to componentToAddFactor
    k = componentToAddFactor;
    if k>-1
        ks = num2str(k);
        %find the data subset for which the posterior is greatest for this component, use tLog
        eval(['[resids,dumLambda,newPsi]=softresidual(data,Lambda',ks,',Psi', ks,',tLogN(:,k),Mu',ks,');']);%newPsi is dummy...
        newFactor = dumLambda(:,size(dumLambda,2));
        %prepare new dumLambda by adding the new factor to Lambda after Lambda_k
        t = sum(numFactors(1,1:k));
        dumLambda = [Lambda(:,1:t) newFactor Lambda(:,t+1:sum(numFactors))];
        dumNumFactorsFactor = numFactors;
        dumNumFactorsFactor(k) = dumNumFactorsFactor(k) + 1;
        %Psi correction
        dumPsi = Psi;
        dumPsi(:,k) = dumPsi(:,k) - newFactor.^2;
        %threshold it to eliminate negative values
        dumPsi = ((dumPsi>=0.0001).*dumPsi) + ((dumPsi<0.0001)*0.0001);

        %let it converge
        trials = 0;
        found = 0;
        while (trials<maxTrial)
            try
                eval(['[dumLambdaFactor,dumPsiFactor,dumMuFactor,dumPiFactor,dumLogsFactor] = cont_mfa2(data,Mu,dumLambda,dumPsi,Pi,dumNumFactorsFactor,-Inf,maxIter);']);
                totalEM = totalEM + size(dumLogsFactor,2); %add to total number of iterations
                newLogs(2) = loglike3(cv,dumLambdaFactor,dumPsiFactor,dumMuFactor,dumPiFactor,dumNumFactorsFactor);
                %set trials to maxTrial to exit while loop
                trials = maxTrial;
                found = 1; %we had no problems
            catch
                %increment the trial counter
                disp('Warning: cont_mfa2 could not add a factor');
                trials = trials + 1;
            end
        end
        %check whether the trials succeeded
        if ~found
            newLogs(2) = -Inf;
        end
    else
        newLogs(2) = -Inf;
    end

    %select the best likelihood
    split = newLogs(1)>newLogs(2); %1 if split, 0 if add factor
    if split
        chosenComponent = componentToSplit;
    else
        chosenComponent = componentToAddFactor;
    end
    CvLog = max(newLogs);
    if CvLog>oldCvLog
        %record what we are doing
        progress(progressCount,2)= chosenComponent;
        progress(progressCount,3)= split;

        %execute the increment
        k = chosenComponent;
        ks = num2str(chosenComponent);

        if split
            %split the chosen component
            disp('Splitting a component');
            Lambda = dumLambdaSplit;
            Psi = dumPsiSplit;
            Mu = dumMuSplit;
            Pi = dumPiSplit;
            dumLogs = dumLogsSplit;
            numFactors = dumNumFactorsSplit;
            numMeans = numMeans + 1;%update the number of Means
        else
            %add a new factor to the new component
            disp('Adding a factor');
            numFactors(k) = numFactors(k) + 1;
            Lambda = dumLambdaFactor;
            Psi = dumPsiFactor;
            Mu = dumMuFactor;
            Pi = dumPiFactor;
            dumLogs = dumLogsFactor;
            factorCount = factorCount+numFactors(k);
        end %if split
        progress(progressCount,1) = size(dumLogs,2);

        %unravel parameters
        factorCount = 0;
        %correct too small values of Lambda and Psi

        patch1 = (abs(Lambda)<lambdaSmall).* (Lambda<0); %negative small values
        patch2 = (abs(Lambda)<lambdaSmall).* (Lambda>=0); %positive small values
        Lambda = Lambda + patch1*(-lambdaSmall) + patch2*lambdaSmall;
        Psi = Psi + (abs(Psi)<0.0001)*0.0001; %and that's a patch...
        for k=1:numMeans
            ks = num2str(k);
            eval(['Lambda',ks,'=Lambda(:,factorCount+1:factorCount+numFactors(k));']);
            factorCount = factorCount + numFactors(k);
            eval(['Mu',ks,'=Mu(:,k);']);
            eval(['Psi',ks,'=Psi(:,k);']);
            %calculate tLog, weighted (but not normalized to sum up to 1) posteriors
            eval(['tLog(:,k) = Pi(k)*exp(loglike2(data,Lambda',ks,',Psi',ks,',Mu',ks,'));']);
        end

        %normalize tLog
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
                smallLogIndices = find((sum(tLogN')<0.1));
                for i=smallLogIndices
                    [dummy maxEl] = max(tLog(i,:));
                    tLogN(i,maxEl)=1;
                    tLogN(i,:) = tLogN(i,:) / sum(tLogN(i,:));
                end
            end
        end %normalize tLog

        logsTra = [logsTra dumLogs];
        logsCv = [logsCv CvLog];
        progressCount = progressCount+1;
    end %if cvlog>oldlog

    %disp(progress); %if you want to track factor and component addition
end %while

count = 0;
for i = 1:size(numFactors,2)
    mixture{i}.Lambda = Lambda(:,count+1:count+numFactors(i));
    mixture{i}.Mu= Mu(:,i);
    mixture{i}.Pi = Pi(i);
    mixture{i}.Psi = Psi(:,i);
    mixture{i}.numFactors = numFactors(i);
    count = count + numFactors(i);
end

trainingHistory.logsTra = logsTra;
trainingHistory.logsCv = logsCv;
trainingHistory.totalEM = totalEM;
trainingHistory.progress = progress;

