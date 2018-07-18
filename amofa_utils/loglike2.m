function a =loglike2(data,lambda,psi,mu)
%loglike2(data,lambda,psi,mu) returns the log likelihood of the factor analysis model for the given data set ->returns a numSamples x 1 matrix!
% data   numSamples x numDims data
% lambda numDims x numFactors factor loading matrix
% psi    numDims x 1          the diagonal of the noise covariance matrix
% mu     numDims x 1          the mean vector of the noise
%
% 24 september 2003 - albert ali salah

numSamples = size(data,1);
numDims = size(data,2);
numFactors = size(lambda,2);

%find covariance
Cov = lambda*lambda' + diag(psi);
eps = 0.0005;
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

lTerm = -0.5*numDims*log(2*pi)-0.5*(log(dC)-numDims*log(regTerm));
%like = zeros(numSamples,1);

centered=data'-repmat(mu,1,numSamples);
%for i=1:numSamples
%  x = data(i,:)';
    %HK: faster computation via matrix operations
    if numDims~=1
        like = lTerm-.5*(sum(centered.*(Cov\centered),1))';  
    else
        like = lTerm-.5*(Cov\(centered.^2))';
    end
%end

a = like;

end