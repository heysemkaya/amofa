% function [Lambda,Psi,Mu,Pi,logs]=cont_mfa(data,iMu,iLambda,iPsi,iPi,numFactors,iLik,maxIter,tol);
% 
% Maximum Likelihood Mixture of Factor Analysis using EM continued; uses given initial parameters
%
% data - data matrix
% iMu - initial means (numDims x numMeans)
% iLambda - initial factor loadings: numDims x (sum(numFactors))
% iPsi - initial uniquenesses matrix (numDims x numMeans) <Psi are alllowed to be different>
% iPi - initial component priors (numMeans x 1)
% numFactors - number of factors for each component (1 x numMeans)
% iLik - initial log-likelihood
% maxIter - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% Lambda - factor loadings 
% Psi - diagonal uniqueness matrices
% Mu - mean vectors
% Pi - priors
% logs - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or maxIter steps of EM 
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

function [Lambda,Psi,Mu,Pi,numFactors,logs,minDL]=cont_mfa_fj(data,iMu,iLambda,iPsi,iPi,inumFactors,iDL,maxIter,tol)

if nargin<9 
    tol=0.00001; 
end;
if nargin<8   
    maxIter=100; 
end;

gamma=0.0001;
Phmin=0.0001; 
lambdaSmall = 0.001;
factAnnihThreshold=0.1;
[numSamples, numDims]=size(data);
numMeans = size(iMu,2);
tiny=exp(-700);

%initialize Lambda, concatenation of all mofa{k}.Lambda
Lambda= iLambda;
%initialize Psi
Psi= iPsi;
%initialize component priors 
Pi=iPi;
%initialize means 
Mu=iMu;
minDL=iDL;

oldDL = iDL;
dl=iDL;
numFactors=inumFactors;
mofa=cell(numMeans,1);
params=zeros(numMeans,1);
  
logs=[];

H=zeros(numSamples,numMeans); 	% E(w|x) 
factorCount=0;
for k=1:numMeans
 % ks = num2str(k);
  
  mofa{k}.numFactors = numFactors(k);
  mofa{k}.EZ = zeros(numSamples,mofa{k}.numFactors);
  mofa{k}.EZZ = zeros(numSamples,mofa{k}.numFactors);
  mofa{k}.Psi = Psi(:,k);
  mofa{k}.Lambda=Lambda(:,factorCount+1:factorCount+mofa{k}.numFactors);
  mofa{k}.Pi = Pi(k);
  factorCount= factorCount + mofa{k}.numFactors;
  
  %eval(['EZ_',ks,'=zeros(numSamples,mofa{k}.numFactors);']);
  %eval(['EZZ_',ks,'=zeros(mofa{k}.numFactors,mofa{k}.numFactors);']);
end  
XX=zeros(numDims*numMeans,numDims);
s=zeros(numMeans,1);
const=(2*pi)^(-numDims/2);
%%%%%%%%%%%%%%%%%%%%

modl=cell(maxIter,1);
for iter=1:maxIter;

  %%%% E Step %%%%
  factorCount = 0;
  logversions = zeros(1,numMeans); %for the cases we need to estimate H by its log
  flagFactAn=0;
  params=zeros(numMeans,1);
  for k=1:numMeans
    ks = num2str(k);
    I=eye(mofa{k}.numFactors);
    %Psi_ is the inverse of Psi
    Psi_ = diag(1./mofa{k}.Psi);
    %the current component Lambda is mofa{k}.Lambda
    %mofa{k}.Lambda=Lambda(:,factorCount+1:factorCount + mofa{k}.numFactors);
    
    LP=Psi_* mofa{k}.Lambda; %mofa{k}.Lambda;
    auxM1=I+mofa{k}.Lambda'*LP;
    MM=Psi_-LP/auxM1*LP'; %singularity alert> orig : Psi_-LP*
    % determinant regularization patch
    regTerm = 1;    
    detM = det(MM);
    while (detM==0) 
     regTerm = regTerm * 2;
     detM=det(MM*regTerm);
    end
    while (detM==Inf) 
     regTerm = regTerm / 2;
     detM=det(MM*regTerm);
    end
    
    dM = sqrt(detM)*sqrt(regTerm^-numDims);
    
    log_dM = (log(detM)-numDims*log(regTerm))/2;
    Xk=(data-ones(numSamples,1)*Mu(:,k)'); 
    XM=Xk*MM;
    % H are the weighted (by priors Pi) likelihoods
    H(:,k)=log(const*Pi(k))+log_dM-(0.5*rsum(XM.*Xk));
    % we have estimated H(:,k) by its log!!
    logversions(k) = 1;

    %eval(['EZ_',ks,'=XM*mofa{k}.Lambda;']);
    
    mofa{k}.EZ = XM*mofa{k}.Lambda;
    auxNormEZ= sum(mofa{k}.EZ.^2)/numDims;
    indices=auxNormEZ<=factAnnihThreshold;
    %sum(logversions==0) &&
    if ( sum(indices)>0 && mofa{k}.numFactors-sum(indices)>=1) % 
        % if logversions then better to annihilate the whole components
        % later annihilate factor loading vector if the factors are weak
        
        flagFactAn=1;
        disp(['annihilating weak factors for component ' ks]);
        factorsAnnihilated=mofa{k}.Lambda(:,indices);
        %Psi correction
        diagContribution=diag(factorsAnnihilated*factorsAnnihilated');
        mofa{k}.Psi=mofa{k}.Psi+diagContribution;
        %remove corresponding factor loading columns
        mofa{k}.Lambda(:,indices)=[];
        %update numFactors
        numFactors(k)=numFactors(k)-sum(indices);
        mofa{k}.numFactors=mofa{k}.numFactors-sum(indices);
        if (mofa{k}.numFactors>0)
           auxLambda=mofa{k}.Lambda;
            I=eye(mofa{k}.numFactors);
        else
            %if all components are annihilated then specify an artificial
            %factor with little contribution
            auxLambda=ones(numDims,1)*eps;
            I=eps;
        end
        %re-compute expectations
        %Psi_ is the inverse of Psi
        Psi_ = diag(1./mofa{k}.Psi);
        LP=Psi_* auxLambda; 
        auxM1=I+auxLambda'*LP;
        MM=Psi_-LP/auxM1*LP';%singularity alert> orig : Psi_-LP*
        XM=Xk*MM;
        H(:,k)=log(const*Pi(k))+log_dM-(0.5*rsum(XM.*Xk));
        mofa{k}.EZ = XM*auxLambda;
    end
    params(k)=numDims*(mofa{k}.numFactors+2)+1;%+2 for Psi and Means
    
  end
%   if sum(logversions)>0
%       %we have a number of H estimated via log values
%       disp('H contains too small entries. logversions!');
%   end
%   
  
  
  oldDL=dl;
  %check if global Lambda is to be updated
  if (flagFactAn)
      factorCount=0;
      Lambda=zeros(numDims,sum(numFactors));
      for k=1:numMeans
         Lambda(:,factorCount+1:factorCount + mofa{k}.numFactors)=mofa{k}.Lambda; 
         factorCount = factorCount + mofa{k}.numFactors;
      end
      
  end
    
   % H correction for numerical underflow
  Hsum=logsumexp(H,2);
  Hzero=(Hsum==0);  %the entries of H with zero sum
  Nz=sum(Hzero);  %the number of such entries
  H(Hzero,:)=log(tiny/numMeans)*ones(Nz,numMeans);  %replace those with tiny, equal probability
  Hsum(Hzero)=log(tiny)*ones(Nz,1);

  H= H - repmat(Hsum,1,numMeans); %rdiv(H,Hsum); 				
  s_log=logsumexp(H,1); % use logsumexp as intermediate step to eliminate NaN and Inf cond.
  s=exp(s_log);
  s=s+(s==0)*tiny;
  s2=sum(s)+tiny;
  
  
  factorCount = 0;
  for k=1:numMeans  
    ks = num2str(k);
    I = eye(mofa{k}.numFactors);
    Psi_ = diag(1./mofa{k}.Psi);  
    kD=(k-1)*numDims+1:k*numDims;    
    %the current component Lambda is mofa{k}.Lambda
    %mofa{k}.Lambda=Lambda(:,factorCount+1:factorCount + mofa{k}.numFactors);
    factorCount = factorCount + mofa{k}.numFactors;
    
    LP=Psi_*mofa{k}.Lambda;
    auxM1=I+mofa{k}.Lambda'*LP;
     MM=Psi_-LP/auxM1*LP'; %singularity alert > 
    %MM=Psi_-LP*inv(I+mofa{k}.Lambda'*LP)*LP'; %singularity alert
    Xk=(data-ones(numSamples,1)*Mu(:,k)'); 
    XX(kD,:)=rprod(Xk,exp(H(:,k)))'*Xk/s(k); 
    beta=mofa{k}.Lambda'*MM;
    %eval(['EZZ_',ks,'=I-beta*mofa{k}.Lambda +beta*XX(kD,:)*beta'';']); 
    mofa{k}.EZZ= I-beta*mofa{k}.Lambda +beta*XX(kD,:)*beta';
  end;
  
  
  
  
  %lik=sum(log(Hsum+(Hsum==0)*exp(-700)));
  lik = loglike3(data,Lambda,Psi,Mu,Pi,numFactors,1000);
 
  modl{iter}.Lambda=Lambda;
  modl{iter}.Psi=Psi;
  modl{iter}.Mu=Mu;
  modl{iter}.Pi=Pi;
  modl{iter}.lik=lik;
  modl{iter}.numSamples=numSamples;
  modl{iter}.numFactors=numFactors;
  
  
  [dl,~, numparams]=getDescLen(modl{iter},lik,1);
  modl{iter}.dl=dl;
  modl{iter}.numparams=numparams;
  %HK: Log models so as to select the one with maximal loglik upon
  %degeneraiton

 
  
  

  
  %%%% log dl %%%%

  logs=[logs dl];
  %fprintf('cycle %g   \tlog likelihood %g ',iter,lik);
  
  if (iter<=2)
    dlbase=dl;
  elseif (dl>oldDL)
    %HK !!! violation !!!
  elseif abs((dl-oldDL)/(oldDL))<tol 
  break;
  end;

 
  %%%% M Step %%%%
  
  % means and covariance structure
  
  Psi=zeros(numDims,numMeans);
  factorCount = 0;
  sumPsi = zeros(numDims,1);
  
  for k=1:numMeans    
    ks = num2str(k);
    kN=(k-1)*numSamples+1:k*numSamples;

    T0=rprod(data,exp(H(:,k))); %T0 is h_ij*x_i in Zubin's paper
    T1=T0'*[mofa{k}.EZ ones(numSamples,1)];
    %eval(['T1=T0''*[EZ_',ks,' ones(numSamples,1)]; ']);
    %T1 is the first term of eq.15 in Zubin's paper: sum(h_ij*x_i*E(z|x))
    XH=mofa{k}.EZ'*exp(H(:,k));   %orig : eval(['XH=EZ_',ks,'''*H(:,k);']);
    T2= [s(k)*mofa{k}.EZZ XH; XH' s(k)]; %inv([s(k)*mofa{k}.EZZ XH; XH' s(k)]);
    % T2 is inverse of sum of EZZ, but h_ij is taken within the matrix. 
    %for example s(k) is the entry shown with 1 in the EZZ formula, but here it is sum(h_ij) over j only
    
    T3=T1/T2; %T3 is the new augmented Lambda in Zubin's paper
    if sum(sum(isnan(T3)))>0
        disp('T3 has NaN entries... Operation rolled back');
        Lambda= iLambda;
        %initialize Psi
        Psi= iPsi;
        %initialize component priors 
        Pi=iPi;
        %initialize means 
        Mu=iMu;
        minDL=iDL;
        return;
    end
    mofa{k}.Lambda=T3(:,1:mofa{k}.numFactors);
    Lambda(:,factorCount+1:factorCount+mofa{k}.numFactors)=mofa{k}.Lambda;
    factorCount=factorCount+mofa{k}.numFactors;
    Mu(:,k)=T3(:,mofa{k}.numFactors+1);
    
    %T4=diag(T0'*data-T3*T1')/s2;
    %Psi=Psi + T4.*(T4>0);
    mofa{k}.Psi = diag(T0'*data-T3*T1')/s(k);
    %now, sum( Pi(k)*mofa{k}.Psi) will give the good old joint psi...
    sumPsi = sumPsi + Pi(k)*mofa{k}.Psi;
    Psi(:,k)=mofa{k}.Psi;
  end;
  
  
  % numerical stability corrections
  for k=1:numMeans
      patch1 = (abs(mofa{k}.Lambda)<lambdaSmall).* (mofa{k}.Lambda<0); %negative small values
      patch2 = (abs(mofa{k}.Lambda)<lambdaSmall).* (mofa{k}.Lambda>=0); %positive small values
      mofa{k}.Lambda = mofa{k}.Lambda + patch1*(-lambdaSmall) + patch2*lambdaSmall;
      mofa{k}.Psi=(1-gamma)*mofa{k}.Psi+gamma*sumPsi; 
      % if above line contributes the next line is not needed
      mofa{k}.Psi=mofa{k}.Psi.*(mofa{k}.Psi>Phmin)+(mofa{k}.Psi<=Phmin)*Phmin;
  end
  
  Psi = (1-gamma)*Psi+gamma*repmat(sumPsi,1,numMeans); % HK > why lose local manifold info?
  
  Psi=Psi.*(Psi>Phmin)+(Psi<=Phmin)*Phmin; % to avoid zero variances
  
  patch1 = (abs(Lambda)<lambdaSmall).* (Lambda<0); %negative small values
  patch2 = (abs(Lambda)<lambdaSmall).* (Lambda>=0); %positive small values
  Lambda = Lambda + patch1*(-lambdaSmall) + patch2*lambdaSmall;
  
  % priors
  Pi=s'/s2;
  softCount=Pi*numSamples;

  [weakestLink] =getWeakestLink(Pi,numSamples,numDims,numFactors,-1);
  flagCompAnnih=0;
  
      % kill illegitimate components from the weakest on but not the last
      % component
      while weakestLink>0 && numMeans>1
          %if (compInd2Annih(k)==1) || (sum(isnan(mofa{k}.Psi))>0)|| (sum(sum(isnan(mofa{k}.Lambda))) >0)
              k=weakestLink;
              flagCompAnnih=1;
              mofa(k,:)=[];
              disp([ 'annihilating component ' num2str(k)]);
              %this component poses problems, eliminate it...

              params(k)=[];
              Mu(:,k)=[];
              H(:,k)=[];
              numFactors(k)=[];
              Pi(k) =[]; %Pi([1:k-1 k+1:numMeans]);
              Pi = Pi / sum(Pi);
              %update numMeans
              numMeans = numMeans - 1;  
              
              [weakestLink] =getWeakestLink(Pi,numSamples,numDims,numFactors,-1);
%           else
%               k = k + 1;
%           end
          %update factorCount
          
      end   %while
    
  factorCount = 0;
  if (flagCompAnnih)
      Lambda=zeros(numDims,sum(numFactors));
      Psi=zeros(numDims,numMeans);
      for k=1:numMeans
          Lambda(:,factorCount+1:factorCount+mofa{k}.numFactors)=mofa{k}.Lambda;
          Psi(:,k)=mofa{k}.Psi;
          factorCount = factorCount+mofa{k}.numFactors;
      end
      
  end
  
  if (sum(sum(isnan(Psi)))>0)||(sum(sum(isnan(Lambda)))>0)
      disp('Seems we have a problem :)');
  end %if problem...  
  
 
end; %iter
 % when density is >1 in small regions, ML becomes negative!
  [mnLog, minind]=min(logs);
  Lambda= modl{minind,1}.Lambda;
  Psi=modl{minind,1}.Psi;
  Mu=modl{minind,1}.Mu;
  Pi=modl{minind,1}.Pi;
  numFactors=modl{minind,1}.numFactors;
  minDL=modl{minind,1}.dl;

      
end