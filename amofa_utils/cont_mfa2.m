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

function [Lambda,Psi,Mu,Pi,logs,maxLoglik]=cont_mfa2(data,iMu,iLambda,iPsi,iPi,numFactors,iLik,maxIter,tol)

if nargin<9 
    tol=0.0001; 
end;
if nargin<8   
    maxIter=100; 
end;

gamma=0.0001;
Phmin=0.0001; 
lambdaSmall = 0.01;

[numSamples, numDims]=size(data);
numMeans = size(iMu,2);
tiny=exp(-700);

%initialize Lambda, concatenation of all Lambda_k
Lambda= iLambda;
%initialize Psi
Psi= iPsi;
%initialize component priors 
Pi=iPi;
%initialize means 
Mu=iMu;
maxLoglik=iLik;

oldlik = iLik;
lik=iLik;


  
logs=[];

H=zeros(numSamples,numMeans); 	% E(w|x) 

for k=1:numMeans
  ks = num2str(k);
  eval(['EZ_',ks,'=zeros(numSamples,numFactors(k));']);
  eval(['EZZ_',ks,'=zeros(numFactors(k),numFactors(k));']);
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
  for k=1:numMeans
    ks = num2str(k);
    I=eye(numFactors(k));
    %Psi_ is the inverse of Psi
    Psi_ = diag(1./Psi(:,k));
    %the current component Lambda is Lambda_k
    Lambda_k=Lambda(:,factorCount+1:factorCount + numFactors(k));
    factorCount = factorCount + numFactors(k);
    LP=Psi_*Lambda_k;
    MM=Psi_-LP*inv(I+Lambda_k'*LP)*LP'; %singularity alert
    
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
    
    if dM == Inf 
      %disp('dM is infinite in cont_mfa');
      log_dM = (log(detM)-numDims*log(regTerm))/2;
      Xk=(data-ones(numSamples,1)*Mu(:,k)'); 
      XM=Xk*MM;
      % H are the weighted (by priors Pi) likelihoods
      H(:,k)=log(const*Pi(k))+log_dM-(0.5*rsum(XM.*Xk));
      % we have estimated H(:,k) by its log!!
      logversions(k) = 1;
    else
      %Xk is data - mean_k
      Xk=(data-ones(numSamples,1)*Mu(:,k)'); 
      XM=Xk*MM;
      % H are the weighted (by priors Pi) likelihoods
      H(:,k)=const*Pi(k)*dM*exp(-0.5*rsum(XM.*Xk)); 	
    end

    eval(['EZ_',ks,'=XM*Lambda_k;']);
  end
  if sum(logversions)>0
      %we have a number of H estimated via log values
      %error('H contains too small entries');
  end
  
  oldlik=lik;
  %lik=sum(log(Hsum+(Hsum==0)*exp(-700)));
  lik = loglike3(data,Lambda,Psi,Mu,Pi,numFactors,1000);
  %HK: Log models so as to select the one with maximal loglik upon
  %degeneraiton
  modl{iter}.Lambda=Lambda;
  modl{iter}.Psi=Psi;
  modl{iter}.Mu=Mu;
  modl{iter}.Pi=Pi;
  modl{iter}.lik=lik;
  
  Hsum=rsum(H);
  Hzero=(Hsum==0);  %the entries of H with zero sum
  Nz=sum(Hzero);  %the number of such entries
  H(Hzero,:)=tiny*ones(Nz,numMeans)/numMeans;  %replace those with tiny, equal probability
  Hsum(Hzero)=tiny*ones(Nz,1);
  
  H=rdiv(H,Hsum); 				
  s=csum(H);
  s=s+(s==0)*tiny;
  s2=sum(s)+tiny;
  
  factorCount = 0;
  for k=1:numMeans  
    ks = num2str(k);
    I = eye(numFactors(k));
    Psi_ = diag(1./Psi(:,k));  
    kD=(k-1)*numDims+1:k*numDims;    
    %the current component Lambda is Lambda_k
    Lambda_k=Lambda(:,factorCount+1:factorCount + numFactors(k));
    factorCount = factorCount + numFactors(k);
    
    LP=Psi_*Lambda_k;
    MM=Psi_-LP*inv(I+Lambda_k'*LP)*LP'; %singularity alert
    Xk=(data-ones(numSamples,1)*Mu(:,k)'); 
    XX(kD,:)=rprod(Xk,H(:,k))'*Xk/s(k); 
    beta=Lambda_k'*MM;
    eval(['EZZ_',ks,'=I-beta*Lambda_k +beta*XX(kD,:)*beta'';']); 
  
  end;
  
  %%%% log likelihood %%%%

  logs=[logs lik];
  %fprintf('cycle %g   \tlog likelihood %g ',iter,lik);
  
  if (iter<=2)
    likbase=lik;
  elseif (lik<oldlik)
    %HK !!! violation !!!
    %break;
    %fprintf(' violation');
    %elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)|~finite(lik)) 
  elseif abs((lik-oldlik)/(oldlik))<tol 
     
      %lik=models{k}.lik;
  break;
  end;

 % fprintf('\n');
  
  %%%% M Step %%%%
  
  % means and covariance structure
  
  Psi=zeros(numDims,numMeans);
  factorCount = 0;
  sumPsi = zeros(numDims,1);
  
  for k=1:numMeans    
    ks = num2str(k);
    kN=(k-1)*numSamples+1:k*numSamples;

    T0=rprod(data,H(:,k)); %T0 is h_ij*x_i in Zubin's paper
    
    eval(['T1=T0''*[EZ_',ks,' ones(numSamples,1)]; ']);
    %T1 is the first term of eq.15 in Zubin's paper: sum(h_ij*x_i*E(z|x))
    eval(['XH=EZ_',ks,'''*H(:,k);']);
    eval(['T2=inv([s(k)*EZZ_',ks,' XH; XH'' s(k)]); ']);%singularity alert
    % T2 is inverse of sum of EZZ, but h_ij is taken within the matrix. 
    %for example s(k) is the entry shown with 1 in the EZZ formula, but here it is sum(h_ij) over j only
    
    T3=T1*T2; %T3 is the new augmented Lambda in Zubin's paper
    if sum(sum(isnan(T3)))>0
        disp('T3 has NaN entries... Operation rolled back');
        Lambda= iLambda;
        %initialize Psi
        Psi= iPsi;
        %initialize component priors 
        Pi=iPi;
        %initialize means 
        Mu=iMu;
        maxLoglik=iLik;
        return;
    end
    
    Lambda(:,factorCount+1:factorCount+numFactors(k))=T3(:,1:numFactors(k));
    factorCount=factorCount+numFactors(k);
    Mu(:,k)=T3(:,numFactors(k)+1);
    
    %T4=diag(T0'*data-T3*T1')/s2;
    %Psi=Psi + T4.*(T4>0);
    Psi(:,k) = diag(T0'*data-T3*T1')/s(k);
    %now, sum( Pi(k)*Psi(:,k)) will give the good old joint psi...
    sumPsi = sumPsi + Pi(k)*Psi(:,k);
  end;
  Psi = (1-gamma)*Psi+gamma*repmat(sumPsi,1,numMeans);% why lose local manifold info?
  
  Psi=Psi.*(Psi>Phmin)+(Psi<=Phmin)*Phmin; % to avoid zero variances
  patch1 = (abs(Lambda)<lambdaSmall).* (Lambda<0); %negative small values
  patch2 = (abs(Lambda)<lambdaSmall).* (Lambda>=0); %positive small values
  Lambda = Lambda + patch1*(-lambdaSmall) + patch2*lambdaSmall;
  
  % priors
  Pi=s'/s2;
  if (sum(sum(isnan(Psi)))>0)|(sum(sum(isnan(Lambda)))>0)
      disp('problem...');
      k=1;
      factorCount = 1;
      while k<=numMeans
          if (sum(isnan(Psi(:,k)))>0)| ( sum(sum(isnan(Lambda(:,factorCount:numFactors(k))))) >0) 
              %this component poses problems, eliminate it...
              %update Psi
              Psi(:,k:numMeans-1) = Psi(:,k+1:numMeans);
              Psi = Psi(:,1:numMeans-1);
              %update Mu
              Mu(:,k:numMeans-1) = Mu(:,k+1:numMeans);
              Mu = Mu(:,1:numMeans-1);
              %update lambda and numFactors
              if k==1
                  Lambda = mslide(Lambda,sum(numFactors(1:k-1))+1,numFactors(k));
                  numFactors = numFactors(2:size(numFactors,2));
              else
                  Lambda = mslide(Lambda,sum(numFactors(1:k-1))+1,numFactors(k));
                  numFactors(k:size(numFactors,2)-1) = numFactors(k+1:size(numFactors,2));
                  numFactors = numFactors(1:size(numFactors,2)-1);
              end
              %update Pi
              Pi = Pi([1:k-1 k+1:numMeans]);
              Pi = Pi / sum(Pi);
              %update numMeans
              numMeans = numMeans - 1;              
          else
              k = k + 1;
          end
          %update factorCount
          factorCount = sum(numFactors(1:k-1))+1;
      end   %while
  end %if problem...  
 
end; %iter

  [mxLog, maxind]=max(logs);
  Lambda= modl{maxind,1}.Lambda;
  Psi=modl{maxind,1}.Psi;
  Mu=modl{maxind,1}.Mu;
  Pi=modl{maxind,1}.Pi;
  maxLoglik=mxLog;
      
end