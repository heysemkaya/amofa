% function [L,Ph,LL]=ffa(X,numFactors,cyc,tol);
% 
% Fast Maximum Likelihood Factor Analysis using EM
%
% X - data matrix
% numFactors - number of factors
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% L - factor loadings 
% Ph - diagonal uniquenesses matrix
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM 
%
% This code is originally written by Z. Ghahramani, my additions and
% changes are indicated with AAS - Albert Ali Salah

function [L,Ph,LL]=ffa(X,numFactors,cyc,tol);

if nargin<4  tol=0.00001; end;
if nargin<3  cyc=100; end;

initWithPCA=1;

N=length(X(:,1));
numDims=length(X(1,:));
tiny=exp(-700);

%subtract mean from the data
X=X-(ones(N,1)*mean(X));

%XX is X^2 divided by the number of samples
XX=X'*X/N; %XX is approximately cov(X), but divided by N, instead of N-1. thus XX*N/N-1 = cX

% diagXX is the diagonal elements of XX, [sum(X_{i1}^2) sum(X_{i2}^2) .. sum(X_{im}^2)] /N
% with i from 1 to numSamples and m from 1 to numDims
% hence it is (average of squares of first elements, average of squares of second elements...)
diagXX=diag(XX);

cX=cov(X); %the covariance matrix
scale=det(cX)^(1/numDims);
if scale == 0 
    scale = 0.001; 
end %AAS - patch
if scale == Inf 
    scale = realmin^-1;
end %HK - patch

scale = real(scale); %AAS - patch
% initial factor loadings are taken from a zero-mean Gaussian random distribution.
% the covariance however, is scaled by a factor, namely (scale/numFactors).
% because cov(X*f) = cov(X)*f*f
%rng('default');
if (initWithPCA==1)
    [H D] = eig(cX);
    %sort the eigenvalues
     dg= diag(D);
    [Y I] = sort(dg,'descend');%now its is descending
    
    indexMatrix = I(1:numFactors);

    L = H(:,indexMatrix).* dg(indexMatrix).^.5;  
else
    L=randn(numDims,numFactors)*sqrt(scale/numFactors);
end
%Ph is the diagonal entries of the data covariance... 
% he is initializing Ph from the data covariance!
Ph=diag(cX);
Ph = Ph + (Ph<0.00001)*0.00001; %no entry is smaller than 0.00001 AAS

I=eye(numFactors);

lik=0; LL=[];

%to use in the calculation, log of 2pi^(-2/d)
const=-numDims/2*log(2*pi);

for i=1:cyc;

  %%%% E Step %%%%
  
  Phd=diag(1./Ph);%Phd is the inverse of phi  
  LP=Phd*L;%inv(Psi)*L
  % by the matrix inversion lemma, MM is the inverse of the covariance matrix of the data; 
  % M = inv(LL'+Psi)
  MM=Phd-LP*inv(I+L'*LP)*LP';
  
  % AAS - patch for Inf or 0 determinants
  regTerm = 1;
  dM=det(MM*regTerm);
  while (dM==0) 
     regTerm = regTerm * 2;
     dM=det(MM*regTerm);
  end
  while (dM==Inf) 
     regTerm = regTerm / 2;
     dM=det(MM*regTerm);
  end
  if dM == 0 disp('Determinant zero in FFA!!'); end
  logdM = N*.5*log(dM) - N*.5*numDims*log(regTerm);
  %beta is the for the linear projection, beta = L' * inv(phi + LL'). 
  % the expected value of z given x is beta*x.
  
  beta=L'*MM;
  XXbeta=XX*beta';
  EZZ=I-beta*L +beta*XXbeta;

  %%%% Compute log likelihood %%%%
  
  oldlik=lik;
  lik=N*const+logdM-0.5*N*sum(diag(MM*XX));
  %fprintf('cycle %i lik %g \n',i,lik);
  LL=[LL lik];
  
  %%%% M Step %%%%

  L=XXbeta*inv(EZZ);
  Ph=diagXX-diag(L*XXbeta');
  Ph = Ph + (Ph<0.0001)*0.0001; %no entry is smaller than 0.00001 AAS
  
  if (i<=2)    
    likbase=lik;
  elseif (lik<oldlik)     
    disp('VIOLATION');
    %break;
  %elseif ((lik-likbase)<(1+tol)*(oldlik-likbase)|~finite(lik))  
  %elseif ((5)>(lik-oldlik)|~finite(lik))  
  elseif abs((lik-oldlik)/(lik))<0.0001
    break;
  end;

end
