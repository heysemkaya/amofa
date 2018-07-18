function a =loglike3(data,lambda,psi,mu,priors,numFactors,regTerm)
%loglike3(data,lambda,psi,mu,priors,numFactors,regTerm) returns the log likelihood of the MFA for the 
% given data set 
% data   numSamples x numDims data
% lambda numDims x (sum(numFactors) factor loading matrix
% psi    numDims x numMeans         the diagonal of the noise covariance matrix
% mu     numDims x numMeans         the mean vector of the noise
% priors     1 x numMeans           the component priors
% numFactors 1 x numMeans           the number of factors for each component
% regTerm    1 x 1                  against zero determinants, it cancels itself out
%
% 2 september 2003 - albert ali salah

numSamples = size(data,1);
numDims = size(data,2);
numMeans = size(mu,2);
sumLike = 0;
logtiny = -700;
if nargin<7 regTerm = 1; end;
%PATCH 
if min(priors)==0
    %find the entries that are zero, set it to some small number, say 0.05
    priors = (priors==0)*0.05 + priors;
    %now the sum is greater than 1, renormalize
    priors = priors / sum(priors);
end

like = zeros(numSamples,numMeans);
factorCount = 0;
for k = 1:numMeans
  lambda_k = lambda(:,factorCount+1:factorCount+numFactors(k));
  factorCount = factorCount+numFactors(k);
  %find covariance
  Cov = lambda_k*lambda_k' + diag(psi(:,k));
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
  detCov = sqrt(dC)*sqrt(regTerm^-numDims);
  
  pd = priors(k)*(detCov)^-1;
  centered=data'-repmat(mu(:,k),1,numSamples);
 % for i=1:numSamples
  %  x = data(i,:)';
  if numDims ~= 1
   like(:,k) = log(pd)+(-.5*sum(centered.*(Cov\centered)))';  
  else
    like(:,k) = log(pd)-.5*(Cov\(centered.^2))';  
  end
    
  %end
end

like=logsumexp(like,2);

like(like<logtiny) = logtiny;



sumLike = sum(like) - numSamples*numDims*.5*log(2*pi);
a = sumLike;
