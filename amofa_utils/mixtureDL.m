% returns the description lenght in terms of MML (ModelSelCri=1) or MDL/BIC 

function [descLength, numparams]=mixtureDL(mixture,N,modelLogLik,ModelSelCri)

    % 2DO: check ismember()
    numMeans = numel(mixture);

    numparams=0;
    Pi=zeros(numMeans,1);

    for i=1:numMeans
      Pi(i)=  mixture{i}.Pi;
      numparams=numparams+numel(mixture{i}.Lambda)+ ...
      numel(mixture{i}.Mu)+numel(mixture{i}.Pi)+size(mixture{i}.Lambda,2);%last 1 is for #factors
    end

    numparams=numparams+numel(mixture{1}.Psi);%common Psi
    nparsover2= numparams/2;

    if (ModelSelCri==1)
         descLength = -modelLogLik + (nparsover2*sum(log(Pi))) + ...
         (nparsover2 + 0.5)*(numMeans)*log(N/12);

    else 
            descLength=-modelLogLik + nparsover2*log(N);
    end
              
end