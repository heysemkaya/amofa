function [Lambda,Psi,Mu,Pi,numFactors,logs,minDL]=softsplit_fj(data,lambda,psi,posteriors,mu)
%[X1,X2,g,loglik]=softsplit(data,lambda,psi,posteriors,mu) returns two means, their mixture priors & loglik.
% the split is 1 eigenvector to each side of the mean in the PCA direction
% data   numSamples x numDims data
% lambda numDims x numFactors factor loading matrix
% psi    numDims x 1          the diagonal of the noise covariance matrix
% mu     numDims x 1          the mean vector of the noise
% posteriors numSamples x 1  the posteriors of belonging to this component
%
% 8.10.2003 albert ali salah
% updates by heysem kaya: local 2-component MoFA fitting, using a weighted
% sum of all eigenvectors to initialze the new components in the local
% fitting
%use KMeans to initialize means
useKMeans=0;
useMoFA=1;
useWtEig=1;
numDims = size(data,2);
if nargin<5  mu = sum(data.*repmat(posteriors,1,numDims)) / sum(posteriors); end; %soft mean

numFactors=size(lambda,2);
loglik1 = loglike3(data,lambda,psi,mu,1,numFactors);

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
eigval=diag(D);
%sort the eigenvalues
[Y I] = sort(eigval,'descend');%now its is descending
%I = flipud(I); % now the indices are descending
%no_zero_eigs=sum(Y>1e-6);
%eigvectorx=H(:,I(1))*Y(1)^0.5;
% here we compute the weighted sum of all eigenvectors
eigvector = sum(H(:,I)*diag(Y.^0.5),2);  
% X1 and X2 lie to both sides of the mean on the weighted principal direction
%X1 = mu+eigvectorx';
%X2 = mu-eigvectorx';

X1 = mu+eigvector';
X2 = mu-eigvector';

if useKMeans
    [M, v]= kMClustering(data,2,[X1;X2]);
    X1=M(1,:);
    X2=M(2,:);
end

maxIter=100; 
Mu=[X1' X2'];
Pi=[0.5;0.5];
Lambda=[lambda lambda];
Psi=[psi psi];

numFactors=[numFactors numFactors];

if useMoFA

    [Lambda,Psi,Mu,Pi,numFactors,logs,minDL]=cont_mfa_fj(data,Mu,Lambda,Psi,Pi,numFactors,Inf,maxIter);
else
    logs=[];
    Loglik = loglike3(data,Lambda,Psi,Mu,Pi,numFactors);
    model.numSamples=numSamples;
    model.Mu =Mu;
    model.Lambda=Lambda;
    model.Pi=Pi;
    model.Psi=Psi;
    model.numFactors=numFactors;
    minDL=getDescLen(model,Loglik,1);
%loglike3(data,lambda,psi,mu,priors,numFactors,regTerm)
end
% loglike3(data,lambda,psi,mu,priors,numFactors,regTerm)
% X1=X1';
% X2=X2';
% 
% %calculate their likelihoods
% t1 = loglike2(data,lambda,psi,X1);
% t2 = loglike2(data,lambda,psi,X2);
% %calculate priors and add likelihoods
% t_compare = t1>t2;
% g(1,1) = sum(t_compare)/size(data,1);
% g(2,1) = 1-g(1,1);
% %eliminate log of zero cases
% t1 = ((t1<=-700)*(-700))+ (t1>-700).*t1;
% t2 = ((t2<=-700)*(-700))+ (t2>-700).*t2;  
% %add the log of priors to loglikelihoods
% t1=t1 + log(g(1,1));
% t2=t2 + log(g(2,1));
% % HK > PATCHed this part for numeric stability with logsumexp()
% T=[t1 t2];
% %loglik = sum(log(  g(1,1)*exp(t1) + g(2,1)*exp(t2)  ));
% loglik = sum(logsumexp(T,2));