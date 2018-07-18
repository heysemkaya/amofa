function [mixture,trainingHistory] =amofa(data)
% [mixture,trainingHistory] =imofa_mml_ext(data,progressLimit)
% mixture has fields Lambda, Psi, Mu, Pi, numFactors
% numFactors indicates number of factors per component
% trainingHistory has fields logsTra, bestdl, descriptionLength, progress,
% totalEM, models (MoFA models after each adaptive step)

% Original IMoFA-L is developed by Albert Ali Salah
% Extension to AMoFA by Heysem Kaya 
% 1-factor analytic solution based on Z. Ghahramani's ffa code

% -----------------------------------------------------------------------
% Copyleft (2014): Heysem Kaya and Albert Ali Salah
%
% This software is distributed under the terms
% of the GNU General Public License Version 3
% 
% Permission to use, copy, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------

ModelSelCri=1;
maxIter = 100; %max. EM iterations per run
maxTrial = 3; %try cont_mfa 3 times before giving up
%continue incremental model evolution till last 'progressLimit' steps fail
% progressLimit = 1 (no bias), stop as soon as improvement stops
progressLimit=1;

noProgress=0;
[numSamples,numDims] = size(data);
lambdaSmall = 0.001;
dd = numDims * (numDims+2); %for Mardia kurtosis metric
myepsilon=1e-5;
maxFactor= numDims*.8;
downsizeOp=1;% option 1 is min Pi, option 2 is max DL to remove
numFactors = 1; %initial number of factors per component
numMeans = 1;
totalEM = 0; %the total number of EM iterations
logc=log2(2.865064);
%Rissanen J., Information and Complexity in Statistical Modeling, Springer, 2007, p.15

%fit a 1-factor, 1-component (factor analysis) model
[Lambda1,Psi1,logsTra] = ffa3(data,numFactors,maxIter);
Mu1=mean(data);
Pi=1;
totalEM = totalEM + size(logsTra,2); %add to total number of iterations

Mu1 = Mu1'; %means are stored as column vectors

count=1; 
%oldCvLog = -Inf;

DataLogLike = loglike3(data,Lambda1,Psi1,Mu1,1,numFactors,1000);
logsCv = DataLogLike; %keeps track of validation loglike during the whole algorithm
tLog = exp(loglike2(data,Lambda1,Psi1,Mu1)); %the total weighted likelihoods for each sample
%normalize tLog
tLogN = ones(numSamples,1);%normalized posteriors
 
%create the initial batch variables
Lambda = Lambda1;
Psi = Psi1;
Mu = Mu1;

model_0.numSamples=numSamples;
model_0.Mu=Mu;
model_0.Lambda=Lambda;
model_0.Pi=Pi;
model_0.Psi=Psi;
model_0.numFactors=numFactors;
[dlength]=getDescLen(model_0,DataLogLike,ModelSelCri);

dl(1) = dlength;

models{1}.Mu=Mu;
models{1}.Lambda=Lambda;
models{1}.Pi=Pi;
models{1}.Psi=Psi;
models{1}.numFactors=numFactors;
models{1}.dl=dlength;
models{1}.logTra=[];
models{1}.action=-1;

progress(1,1) = size(logsTra,2); %1st column keeps track of iterations
progress(1,2) = 0; %which component
progress(1,3) = 0; %what is done
progressCount = 2;

cont_EM=1;%flag to indicate EM continuation
direction=1;%1 -> up, -1 -> down

while   cont_EM %CvLog > oldCvLog
    action=-1;%No action
    factorCount = 0;
    if numMeans == 1
        componentToSplit = 1;
        componentToAddFactor = 1;
        splitThreshold=numDims*(numFactors(1)+2)+2;
        %in case no split is possible don't attempt to do it!
       if (numSamples < splitThreshold)
           componentToSplit=-1;
       end
        
    else
     
        %
        factorCount2 = 0;
        clear factorMetric
        clear splitMetric
        [ll indic]=max(tLogN,[],2);
            
         %
        if (direction<0)
            
            componentToSplit = -1;
            componentToAddFactor = -1;
            %kill the weakest component
           
           % or kill the component whose message cost is maximum
           factorCount3=0;
           comp_DLs=zeros(numMeans,1);
           for comp=1:numMeans
                cluster_k=(indic==comp);
                localdata=data(cluster_k,:);                
                rngLambda=factorCount3+1:factorCount3+numFactors(comp);
                factorCount3 = factorCount3+numFactors(comp);
                comp_model.numSamples=numSamples;
                comp_model.Mu=Mu(:,comp);
                comp_model.Lambda=Lambda(:,rngLambda);
                comp_model.Pi=Pi(comp);
                comp_model.Psi=Psi(:,comp);
                comp_model.numFactors=numFactors(comp);
                %loglik_comp=loglike3(localdata,Lambda(:,rngLambda),Psi(:,comp),Mu(:,comp),Pi(comp),numFactors(comp),1);
                comp_DLs(comp)=Pi(comp)*numSamples-(numDims*(numFactors(comp)+2))/2;  %getDescLen(comp_model,loglik_comp,1);
                
           end
           
           if (downsizeOp==1)
               [minPi, k]=min(Pi);
           else
               [maxDL, k]=min(comp_DLs);
           end
           
            startL=1;
            endL=numFactors(k);
            if k>1
                startL=sum(numFactors(1:k-1))+1;
                endL=sum(numFactors(1:k));
            end
            
            Mu(:,k)=[];
            Psi(:,k)=[];
            Lambda(:,startL:endL)=[];
            Pi(k)=[];
            Pi=Pi./sum(Pi);
            tLogN(:,k)=[];
            tLogN=tLogN./repmat(sum(tLogN,2),1,size(tLogN,2));
            tLog(:,k)=[];
            numFactors(k)=[];
            numMeans=numMeans-1;
            [ll indic]=max(tLogN,[],2);
            
            disp('Downsizing component annihilation');
            iMu=Mu; iLambda=Lambda; iPsi=Psi;iPi=Pi; inumFactors=numFactors;
            % now run locally detoxified factor analyzers in the mixture
            [Lambda,Psi,Mu,Pi,numFactors,logs,minDL]=cont_mfa_fj(data,iMu,iLambda,iPsi,iPi,inumFactors,Inf,10);
            if (numMeans~=size(Mu,2))
                numMeans=size(Mu,2);
            end
             tLog = exp(loglike4(data,Lambda,Psi,Mu,Pi,numFactors));
             tLogN=tLog./repmat(sum(tLog,2),1,size(tLog,2));
             
            [models]=logModel(models,Mu,Lambda,Pi,Psi,numFactors,dumLogs,minDL,action,tLog,tLogN);
            progress(progressCount,1) = size(dumLogs,2);
            dl(progressCount)=minDL;
            progressCount = progressCount+1;
        else
        
            [ll indic]=max(tLogN,[],2);
            for k=1:numMeans
                ks = num2str(k);
                %find the data set that is associated with this component

                [ll indic]=max(tLogN,[],2);
                cluster_k=(indic==k);
                dat=data(cluster_k,:);


                kD=factorCount2+1:factorCount2+numFactors(k);
                factorCount2 = factorCount2+numFactors(k);
                tLambda = Lambda(:,kD);
                tempCov = tLambda*tLambda' +diag(Psi(:,k));

                %Mardia multivariate kurtosis metric
                invCov = inv(tempCov);
                tMu = Mu(:,k)';
                b_sum = 0;
                for i=1:numSamples
                    x = data(i,:) - tMu;%centered data
                    b_sum = b_sum + tLogN(i,k)*(x*invCov*x')^2;
                end
                b_sum = b_sum / sum(tLogN(:,k));
                splitMetric(k) = abs((b_sum - dd) / sqrt(8*dd/sum(tLogN(:,k))));

                % factor metric: the sum of squares of the covariance differences
                factorMetric(k) = sum(sum(((tempCov - cov(dat)).*(1-eye(size(tempCov,1)))).^2));

                %if the size of the component is too small, do not split it
                % here we check if the equal sized newborn components will
                % survive the next step
              
                % component params *2 
                % for numFactors consider the new case case
                % (numfactors(k)+1+2)
                factorThreshold=(numDims*(numFactors(k)+3)+getIntCodeLen(numFactors(k)+1)+logc)/2;
                % for a component to split it needs twice support of its
                % current annihilation threshold: C(k)=T(k)/2 so threshold
                % is C(k) plus the cost for mixture proportions
                splitThreshold=numDims*(numFactors(k)+2)+2; 

                if Pi(k)*numSamples <= splitThreshold
                    disp(['Split not allowed due to insufficient soft support for component ' ks]);
                    splitMetric(k) = -Inf;
                end
                %if the number of factors are too much, do not consider adding new factors
                if numFactors(k)>numDims*.8 || Pi(k)*numSamples <= factorThreshold
                    factorMetric(k) = -Inf;
                end
            end
            %we have looped for all components

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
    end
    
    if componentToAddFactor>0 && numFactors(componentToAddFactor)>=maxFactor
            componentToAddFactor = -1;
    end
    
    %initialize newLogs, to select which action pays off better
    newLogs = zeros(1,2);

    %split the componentToSplit
    k = componentToSplit;
    if k>-1 && direction>0
        ks = num2str(k);
        %first find the split
        [ll indic]=max(tLogN,[],2);
        cluster_k=(indic==k);
        localdata=data(cluster_k,:);
        if numel(localdata)<=1
           disp('No local data for splitting'); 
           componentToSplit=-1;
           newLogs(1) = -Inf;
           newDlength(1)=Inf;
        else
            eval(['[Lambda_l,Psi_l,Mu_l,Pi_l,NumFactors_l,logs_l,minDL_l] = softsplit_fj(localdata,Lambda',ks,',Psi',ks,',tLogN(cluster_k,k),Mu',ks,');']);
            %create a dummy set of batch variables
            oldFactorCount=sum(numFactors);
            factorCount = sum(numFactors(1,1:k-1));
            found = 0;
            reject=0;
            if (minDL_l~=Inf)
                if (size(Mu_l,2)==2)

                    % update Lambda: the local dimensionality of new components
                    % may vary

                    dumLambda = [Lambda(:,1:factorCount) Lambda_l(:,1:NumFactors_l(1)) ...
                        Lambda(:,factorCount+numFactors(k)+1:oldFactorCount) Lambda_l(:,NumFactors_l(1)+1:sum(NumFactors_l))];

                    tmpNumFactors=numFactors;

                    tmpNumFactors(k)=NumFactors_l(1);%local split can adapt the number of factors :)
                    dumNumFactorsSplit = [tmpNumFactors NumFactors_l(2)];
                    %Init components for global EM sweep
                    dumPsi = [Psi Psi_l(:,2)];
                    dumPsi(:,k)=Psi_l(:,1);
                    dumMu = [Mu Mu_l(:,2)];
                    dumMu(:,k) = Mu_l(:,1);
                    dumPi = Pi;
                    prevPi=dumPi(k,1) ;
                    dumPi(k,1) = prevPi *Pi_l(1);  %the prior of the component is halved
                    dumPi(numMeans+1,1) = prevPi*Pi_l(2);

                else
                    % if one of the components are annihilated during
                    % initialization of split then use the new params
                    reject=1;
                    disp(['Split is rejected during initialization for component ' ks]);
                    dumMu(:,k)=Mu_l(:,1);
                    dumNumFactorsSplit(k)=NumFactors_l(1);
                    dumLambda = [Lambda(:,1:factorCount) Lambda_l(:,1:NumFactors_l(1)) ...
                        Lambda(:,factorCount+numFactors(k)+1:oldFactorCount)];

                    dumPsi(:,k)=Psi_l(:,1);
                    % no change in dumPi
                    dumPi = Pi;
                end
                
                if (reject==0)
                    trials = 0;

                    while (trials<maxTrial)
                        try
                            [dumLambdaSplit,dumPsiSplit,dumMuSplit,dumPiSplit,dumNumFactorsSplit,dumLogsSplit,newDlength(1)] =...
                                cont_mfa_fj(data,dumMu,dumLambda,dumPsi,dumPi,dumNumFactorsSplit,Inf,maxIter);
                            totalEM = totalEM + size(dumLogsSplit,2); %add to total number of iterations

                            model_cs.numSamples=numSamples;
                            model_cs.Mu=dumMuSplit;
                            model_cs.Lambda=dumLambdaSplit;
                            model_cs.Pi=dumPiSplit;
                            model_cs.Psi=dumPsiSplit;
                            model_cs.numFactors=dumNumFactorsSplit;
                            model_cs.totalEM=totalEM;

                            trials = maxTrial;
                            found = 1; %we had no problems
                        catch err
                            %increment the trial counter
                            disp(['Warning: cont_mfa2 could not split component ' ks]);
                            trials = trials + 1;
                        end
                    end
                else
                    newDlength(1)=Inf;
                end
                
            end
            %check whether the trials succeeded
            if ~found
                %newLogs(1) = -Inf;
                newDlength(1)=Inf;
            end
        end
    else
        %newLogs(1) = -Inf;
        newDlength(1)=Inf;
    end

    %now try adding a factor to componentToAddFactor
    k = componentToAddFactor;
    if k>-1 && direction>0
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
                [dumLambdaFactor,dumPsiFactor,dumMuFactor,dumPiFactor,dumNumFactorsFactor,dumLogsFactor,newDlength(2)] = ...
                    cont_mfa_fj(data,Mu,dumLambda,dumPsi,Pi,dumNumFactorsFactor,Inf,maxIter);
                totalEM = totalEM + size(dumLogsFactor,2); %add to total number of iterations
                %newLogs(2) = loglike3(data,dumLambdaFactor,dumPsiFactor,dumMuFactor,dumPiFactor,dumNumFactorsFactor);

                model_f.numSamples=numSamples;
                model_f.Mu=dumMuFactor;
                model_f.Lambda=dumLambdaFactor;
                model_f.Pi=dumPiFactor;
                model_f.Psi=dumPsiFactor;
                model_f.numFactors=dumNumFactorsFactor;
                model_f.totalEM=totalEM;
               
                trials = maxTrial;
                found = 1; %we had no problems
            catch err
                %increment the trial counter
                disp(['Warning: cont_mfa2 could not add a factor to component' ks]);
                trials = trials + 1;
            end
        end
        %check whether the trials succeeded
        if ~found
            %newLogs(2) = -Inf;
            newDlength(2)=Inf;
        end
    else
        %newLogs(2) = -Inf;
        newDlength(2)=Inf;
    end
    
    %select the best likelihood
    [minDlenght, action] = min(newDlength); 
    if action==1
        chosenComponent = componentToSplit;
    elseif action==2
        chosenComponent = componentToAddFactor;
%     elseif action==3
%         chosenComponent = componentToAnnihilate;
    end
 
    deltaProgress=abs((minDlenght-dl(progressCount-1))/dl(progressCount-1));
    
    if (minDlenght<dl(progressCount-1)) && (deltaProgress>myepsilon)
        noProgress=0;
    else
        noProgress=noProgress+1;
        if ((direction>0) && (noProgress>=progressLimit || minDlenght==Inf))
                direction=-1;
                noProgress=0;
        elseif (direction < 0 && numMeans==1)
            cont_EM=0;
        end
    end %if cvlog>oldlog
    
        if (minDlenght~=Inf)
            %record what we are doing
            progress(progressCount,2)= chosenComponent;
            progress(progressCount,3)= action;

            %execute the increment
            k = chosenComponent;
            ks = num2str(chosenComponent);

            if action==1 
                %split the chosen component
                disp(['Splitting component ' ks]);
                Lambda =model_cs.Lambda;  % dumLambdaSplit;
                Psi = model_cs.Psi;% dumPsiSplit;
                Mu = model_cs.Mu;% dumMuSplit;
                Pi = model_cs.Pi;% dumPiSplit;
                dumLogs = dumLogsSplit;
                numFactors =model_cs.numFactors; % dumNumFactorsSplit;
                [nm]=size(Mu,2);

                numMeans =nm;  
            elseif action==2
                %add a new factor to the new component
                disp(['Adding a factor to component ' ks]);
                Lambda =model_f.Lambda;  % ;
                Psi = model_f.Psi;% ;
                Mu = model_f.Mu;% ;
                Pi = model_f.Pi;% ;
                numFactors =model_f.numFactors; % dumNumFactorsSplit;
                dumLogs = dumLogsFactor;
                %factorCount = factorCount+numFactors(k);
                [nm]=size(Mu,2);

                numMeans =nm; 
            end %if action...
            [models]=logModel(models,Mu,Lambda,Pi,Psi,numFactors,dumLogs,minDlenght,action,tLog,tLogN);
            progress(progressCount,1) = size(dumLogs,2);
            dl(progressCount)=minDlenght;

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
               
            end
           
            tLog = exp(loglike4(data,Lambda,Psi,Mu,Pi,numFactors));
            %normalize tLog
            if size(tLog,2)==1
                tLogN = ones(numSamples,1);
            else
                tLogN = tLog ./ repmat(sum(tLog,2),1,size(tLog,2));
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
            %logsDL = [logsCv CvLog];
            progressCount = progressCount+1;
        
        end

end %while

count = 0;
[bestDL, bestind]=min(dl);
 bestmodel=models{bestind};

for i = 1:size(bestmodel.numFactors,2)
    mixture{i}.Lambda = bestmodel.Lambda(:,count+1:count+bestmodel.numFactors(i));
    mixture{i}.Mu=  bestmodel.Mu(:,i);
    mixture{i}.Pi =  bestmodel.Pi(i);
    mixture{i}.Psi =  bestmodel.Psi(:,i);
    mixture{i}.numFactors =  bestmodel.numFactors(i);
    count = count +  bestmodel.numFactors(i);
end
trainingHistory.bestdl=bestDL;
trainingHistory.models=models;
trainingHistory.logsTra = logsTra;
trainingHistory.descriptionLength = dl;
trainingHistory.totalEM = totalEM;
trainingHistory.progress = progress;
disp('fini');
% catch err
%  disp(err.identifier);
%     
% end
