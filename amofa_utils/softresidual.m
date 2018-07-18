
function [resids,newLambda,newPsi]=softresidual(data,lambda,psi,posteriors,mu);
%softresidual(data,lambda,psi,posteriors,mu) takes a data set and the parameters of a factor analyser. 
% calculates the expected values of the factors, find the difference of these estimates and the 
% actual data points(called the residuals) and returns the principal component of the residuals 
% this vector will be used as the initial value of the k+1th factor (we are increasing the number of factors)
% and thus concatenated to lambda. psi will be changed accordingly. 
% data        numSamples x numDims data
% lambda      numDims x numFactors factor loading matrix
% psi         numDims x 1          the diagonal of the noise covariance matrix
% posteriors  numSamples x 1       normalized posteriors for each sample and this component
% mu          numDims x 1          the mean vector of the noise
%
% 8.10.2003 albert ali salah
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

numDims = size(data,2);
if nargin<5  mu = sum(data.*repmat(posteriors,1,numDims)) / sum(posteriors); end; %soft mean
if nargin<4  posteriors = ones(size(data,1),1); end;

numSamples = size(data,1);
numFactors = size(lambda,2);

if size(mu,1) == 1
  mu = mu';
end  

% according to Ghahramani & Hinton E[z] is lambda'(psi+lambda*lambda')^-1*x
premult_z = lambda'*inv(diag(psi)+lambda*lambda');

%find the residual data set
for i=1:numSamples
  x = data(i,:)'-mu;
  z = premult_z*x;
  E_x = lambda*z;
  resids(i,:) = (x - E_x)';  
end

%find the soft covariance
resids_mu = sum(resids.*repmat(posteriors,1,numDims)) / sum(posteriors); %mean of resids
c_data= resids - repmat(resids_mu,numSamples,1); %centered residuals
ch_data = c_data .* repmat(posteriors.^.5,1,numDims); %weight each data point by the square root of posterior
softCov = ch_data'*ch_data/(sum(posteriors)-1); %during multiplication, the posterior will be reconstructed

%Patch heysem 19.10.2012
if (sum(sum(isnan(softCov)))>0)
   disp('Warning: softCov contains NaN elems');
   softCov(isnan(softCov))=eps;
   %softCov
end
if (sum(sum(isinf(softCov)))>0)
   disp('Warning: softCov contains inf elems');
   softCov(isinf(softCov))=eps;
   %softCov
end
    
%fit one factor to the residue
[H D] = eig(softCov);
%sort the eigenvalues
[Y I] = sort(diag(D));%default is ascending
I = flipud(I); % now the indices are descending
indexMatrix = I(1);
LambdaR = H(:,indexMatrix)*D(indexMatrix,indexMatrix)^.5;

%combine the two
newLambda = [lambda LambdaR];
newPsi = 0;%dummy