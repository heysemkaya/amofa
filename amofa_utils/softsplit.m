function [X1,X2,g,loglik]=softsplit(data,lambda,psi,posteriors,mu)
%[X1,X2,g,loglik]=softsplit(data,lambda,psi,posteriors,mu) returns two means, their mixture priors & loglik.
% the split is 1 eigenvector to each side of the mean in the PCA direction
% data   numSamples x numDims data
% lambda numDims x numFactors factor loading matrix
% psi    numDims x 1          the diagonal of the noise covariance matrix
% mu     numDims x 1          the mean vector of the noise
% posteriors numSamples x 1  the posteriors of belonging to this component
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

%use KMeans to initialize means
useKMeans=0;

numDims = size(data,2);
if nargin<5  mu = sum(data.*repmat(posteriors,1,numDims)) / sum(posteriors); end; %soft mean

loglik1 = sum(loglike2(data,lambda,psi,mu));


if nargin<4  posteriors = ones(size(data,1),1); end;

if (size(mu,1)>1)
    mu=mu';
end

tiny = exp(-700);%we can use realmin here
g = zeros(2,1);
numSamples = size(data,1);

%find the principal direction by soft covariance
c_data= data - repmat(mu,numSamples,1); %centered data  
ch_data = c_data .* repmat(posteriors.^.5,1,numDims); %weight each data point by the square root of posterior
softCov = ch_data'*ch_data/(sum(posteriors)-1); %during multiplication, the posterior will be reconstructed

%Patch heysem 22.03.2013
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


[H D] = eig(softCov);
%sort the eigenvalues
[Y I] = sort(diag(D),'descend');%now its is descending
%I = flipud(I); % now the indices are descending
indexMatrix = I(1);

eigvector = H(:,indexMatrix)*D(indexMatrix,indexMatrix)^.5;  
% X1 and X2 lie to both sides of the mean on the principal direction
X1 = mu+eigvector';
X2 = mu-eigvector';

if useKMeans
    [M, v]= kMClustering(data,2,[X1;X2]);
    X1=M(1,:);
    X2=M(2,:);
end

X1=X1';
X2=X2';

%calculate their likelihoods
t1 = loglike2(data,lambda,psi,X1);
t2 = loglike2(data,lambda,psi,X2);
%calculate priors and add likelihoods
t_compare = t1>t2;
g(1,1) = sum(t_compare)/size(data,1);
g(2,1) = 1-g(1,1);
%eliminate log of zero cases
t1 = ((t1<=-700)*(-700))+ (t1>-700).*t1;
t2 = ((t2<=-700)*(-700))+ (t2>-700).*t2;  
%add the log of priors to loglikelihoods
t1=t1 + log(g(1,1));
t2=t2 + log(g(2,1));
% HK > PATCHed this part for numeric stability with logsumexp()
T=[t1 t2];
%loglik = sum(log(  g(1,1)*exp(t1) + g(2,1)*exp(t2)  ));
loglik = sum(logsumexp(T,2));