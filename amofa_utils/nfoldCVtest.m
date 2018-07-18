function [folds,PerformParams]=nfoldCVtest(data,datalabels,T,dset,foldCnt)
numC=max(datalabels);%number of classes
%N fold CV
if nargin <5
    foldCnt=10;
end
PerformParams=zeros(foldCnt,2);% performance and parameters are kept
folds=cell(foldCnt,1);
[N, dimens]=size(data);

% dataset name
if nargin <4
    dset='dataset';
end
% AMoFA parameter: number of adaptive iterations without improvement
if nargin <3
    T=1;
end


step=round(N/foldCnt);

for sim=1:foldCnt
    
     sti=(sim-1)*step+1;
     endi=sim*step;
     if sim==foldCnt
         endi=N;
     end
     testIndic=zeros(N,1);
     testIndic(sti:endi)=1;
     
     folds{sim}.tsdata=data(testIndic==1,:);
     folds{sim}.tslabels=datalabels(testIndic==1,:);
     folds{sim}.trdata= data(testIndic==0,:);
     folds{sim}.trdatalabels=datalabels(testIndic==0,:);
    for i=1:numC
        
            folds{sim}.class{i}.trdata=folds{sim}.trdata(folds{sim}.trdatalabels==i,:);
    end 
    
end


for sim=1:foldCnt

   for j=1:numC
       
        [folds{sim}.class{j}.mixture_amofa, folds{sim}.class{j}.trainingHistory_amofa] = ...
        amofa(folds{sim}.class{j}.trdata);
    end

    %evaluate the contributions on all test samples
    
    folds{sim}.Nts=size(folds{sim}.tsdata,1);
    folds{sim}.posts=zeros(folds{sim}.Nts,numC);
    folds{sim}.params=0;
    for k = 1:numC
       %Lambda,Mu,Psi,Pi,numFactors
      [folds{sim}.Lambda,folds{sim}.Mu,folds{sim}.Psi,folds{sim}.Pi,folds{sim}.numFactors] = unpack(folds{sim}.class{k}.mixture_amofa);
      % data,lambda,psi,mu,priors,numFactors,regTerm
      folds{sim}.loglikeVal = loglike4(folds{sim}.tsdata,folds{sim}.Lambda,folds{sim}.Psi,folds{sim}.Mu,folds{sim}.Pi,folds{sim}.numFactors);
      folds{sim}.params=folds{sim}.params+ numel(folds{sim}.Lambda)+numel(folds{sim}.Psi)+numel(folds{sim}.Mu)+numel(folds{sim}.Pi);
      if size(folds{sim}.loglikeVal,2)>1
        %more than one component
        folds{sim}.loglikeVal = logsumexp(folds{sim}.loglikeVal,2);
      end
      folds{sim}.posts(:,k) = folds{sim}.loglikeVal;
    end

    %assign each test sample to the class with highest likelihood
    [dummy, folds{sim}.foundLabels ] = max(folds{sim}.posts,[],2);
    %output classification accuracy
    folds{sim}.accuracy = sum(folds{sim}.tslabels==folds{sim}.foundLabels)*100/max(size(folds{sim}.tslabels));
    disp(['Classification accuracy is ',num2str(folds{sim}.accuracy),' per cent.']);
end


for simi=1:foldCnt
    PerformParams(simi,1)=folds{simi}.accuracy;
    PerformParams(simi,2)=folds{simi}.params;
end
str=strcat(dset,'_results_amofa');
save(str,'PerformParams','folds');

