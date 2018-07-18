function [a, like] =loglike_mix(data,mixture,isGmm)
% [a, like] =loglike_mix(data,mixture,isGmm) returns the log likelihood of the mixture for the 
% given data set : Heysem Kaya
% data   numSamples x numDims data
% mixture which has
% - Lambda numDims x (sum(numFactors) factor loading matrix
% - Psi    numDims x numMeans         the diagonal of the noise covariance matrix
% - Mu     numDims x numMeans         the mean vector of the noise
% - Pi     1 x 1 the mixture proportion
% the 
% isGmm    1 x 1                  specifying whether mixture is a GMM or
% MoFA incase isGmm=1 instead of lambda and psi Cov contains covariance
% info
% % 
% 2 september 2003 - albert ali salah
numMeans = numel(mixture);
numSamples = size(data,1);
numDims = size(data,2);
if nargin<3 
    isGmm=0; 
end
    
priors=zeros(numMeans,1);
for k=1:numMeans
   priors(k)=mixture{k}.Pi;  
end

sumLike = 0;
logtiny = -700;
regTerm = 1; 

%PATCH 
if min(priors)==0
    %find the entries that are zero, set it to some small number, say 0.05
    priors = (priors==0)*0.05 + priors;
    %now the sum is greater than 1, renormalize
    priors = priors / sum(priors);
end

like = zeros(numSamples,numMeans);
%factorCount = 0;
for k = 1:numMeans
    if (~isGmm)
     lambda_k = mixture{k}.Lambda;
     % factorCount = factorCount+numFactors(k);
     %find covariance
     Cov = lambda_k*lambda_k' + diag(mixture{k}.Psi);
    else
     Cov=   mixture{k}.Cov;
    end
  eps = 0.0005;
  
  if sum(sum(isinf(Cov) + isnan(Cov)))>0
    %problem....
    disp('problem...');
  end  
  if rank(Cov)<numDims
    Cov = Cov + eye(numDims)*eps;
  end  

  %invCov = inv(Cov);
  regTerm = 1;
  dC = det(Cov);
  while (dC==0) 
    regTerm = regTerm * 2;
    dC=det(Cov*regTerm);
  end
  while (dC==Inf)
    regTerm = regTerm / 2;
    dC=det(Cov*regTerm);
  end
  %detCov = sqrt(dC)*sqrt(regTerm^-numDims);
  
  %pd = priors(k)*(detCov)^-1;
  coef=log(priors(k)) -(numDims*.5)*log(2*pi) -.5*(log(dC) -numDims*log(regTerm)) ;
  centered=data'-repmat(mixture{k}.Mu,1,numSamples);
 % for i=1:numSamples
  %  x = data(i,:)';
    like(:,k) = coef+(-.5*sum(centered.*(Cov\centered)))';  
  %end
end
like=logsumexp(like,2)-numDims*.5*log(2*pi);
like(like<logtiny) = logtiny;
sumLike = sum(like); %- numSamples*numDims*.5*log(2*pi);
a = sumLike;