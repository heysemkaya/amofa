function a =loglike(data,lambda,psi,mu);
%loglike(data,lambda,psi,mu) returns the log likelihood of the factor analysis model for the given data set ->returns the summation only, for individual points use loglike2
% data   numSamples x numDims data
% lambda numDims x numFactors factor loading matrix
% psi    numDims x 1          the diagonal of the noise covariance matrix
% mu     numDims x 1          the mean vector of the noise
%
% 23 september 2003 - albert ali salah

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
like = 0;
for i=1:numSamples
  x = data(i,:)';
  like = like + (x-mu)'/Cov*(x-mu);  
end

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

a = -.5*(numSamples*(numDims*log(2*pi)+log(dC)-numDims*log(regTerm))+like);  