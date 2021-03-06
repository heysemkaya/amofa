function a =loglike_gmm(data,mu,covM,priors,regTerm)
%loglike4(data,lambda,psi,mu,priors,numFactors,regTerm) returns the log likelihood of the gmm for the 
% given data set - returns numSamples x numComponents matrix that needs to be normalized
% adaptepted from MoFA for GMMs
% data   numSamples x numDims data
% cov numDims x numDims x numMeans covariance matrices
% mu     numDims x numMeans         the mean vector of the noise
% priors     1 x numMeans           the component priors
% regTerm    1 x 1                  against zero determinants, it cancels itself out
%
% 29 jan 2004 - albert ali salah

numSamples = size(data,1);
numDims = size(data,2);
numMeans = size(mu,2);
sumLike = 0;
tiny = exp(-700);
if nargin<5 regTerm = 1; end;

%PATCH --> Heysem: We can simply remove those components from the mixture
if min(priors)==0
    %find the entries that are zero, set it to some small number, say 0.05
    priors = (priors==0)*0.05 + priors;
    %now the sum is greater than 1, renormalize
    priors = priors / sum(priors);
end

    
like = zeros(numSamples,1);
factorCount = 0;
a = zeros(numSamples,numMeans);

for k = 1:numMeans

  %find covariance
  Cov = covM(:,:,k); 
  eps = 0.0001;
  
  if sum(sum(isinf(Cov) + isnan(Cov)))>0
    %problem....
    disp('problem...');
  end  
  if rank(Cov)<numDims
    Cov = Cov + eye(numDims)*eps;
  end  

 % invCov = inv(Cov);
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
  
  coef=log(priors(k)) -(numDims*.5)*log(2*pi) -.5*(log(dC) -numDims*log(regTerm)) ;
  centered=data'-repmat(mu(:,k),1,numSamples);
  %for i=1:numSamples
   % x = data(i,:)';
    a(:,k) = (coef- (.5*sum(centered.*(Cov\centered))))';  
  %end
end

