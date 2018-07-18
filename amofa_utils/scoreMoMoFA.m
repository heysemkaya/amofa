function [accuracy,params,foundLabels]=scoreMoMoFA(class,testdata,groundTruth,numC,priors)

    addpriors=0;
    if (nargin==5)
        logpriors=log(priors);
        if size(logpriors,1)>1
            logpriors=logpriors';
        end
        addpriors=1;
    end
    Nts=size(testdata,1);
    posts=zeros(Nts,numC);
    params=0;
    for i = 1:numC
       %Lambda,Mu,Psi,Pi,numFactors
      [Lambda,Mu,Psi,Pi,numFactors] = unpack(class{i}.mixture);
      % data,lambda,psi,mu,priors,numFactors,regTerm
      loglikeVal = loglike4(testdata,Lambda,Psi,Mu,Pi,numFactors);
      params=params+ numel(Lambda)+numel(Psi)+numel(Mu)+numel(Pi);
      if size(loglikeVal,2)>1
        %more than one component
        loglikeVal = logsumexp(loglikeVal,2);
      end
      posts(:,i) = loglikeVal;
    end
    if (addpriors)
        posts=posts+repmat(logpriors,Nts,1);%has the effect o multiplying with class priors
    end
    %assign each test sample to the class with highest likelihood
  
    [dummy, foundLabels ] = max(posts,[],2);
     accuracy = sum(groundTruth==foundLabels)*100/max(size(groundTruth));   
end