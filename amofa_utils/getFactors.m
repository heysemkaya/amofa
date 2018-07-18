function [Factors,components] = getFactors(data,mixture)
% extract factors of data given mofa
% data: n x d dataset
% mixture: mixture solution outputed by amofa
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
numMeans=numel(mixture);
[N,d]=size(data);
components=cell(numMeans,1);
numFactors=zeros(numMeans,1);
totFactors=0;
Mu=zeros(d,numMeans);

for k=1:numMeans
  totFactors= totFactors+mixture{k}.numFactors;
  numFactors(k)=mixture{k}.numFactors;
  Mu(:,k)=mixture{k}.Mu;
  Pi(k)=mixture{k}.Pi;
end

const=(2*pi)^(-d/2);
Factors=zeros(N,totFactors);
factorCount=0;
  for k=1:numMeans
    ks = num2str(k);
    I=eye(numFactors(k));
    %Psi_ is the inverse of Psi
    Psi_ = diag(1./mixture{k}.Psi);
    %the current component Lambda is Lambda_k
    Lambda_k= mixture{k}.Lambda; %Lambda(:,factorCount+1:factorCount + numFactors(k));
   
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
    
    dM = sqrt(detM)*sqrt(regTerm^-d);
    
    if dM == Inf 
      %disp('dM is infinite in cont_mfa');
      log_dM = (log(detM)-numDims*log(regTerm))/2;
      Xk=(data-ones(N,1)*Mu(:,k)'); 
      XM=Xk*MM;
      % H are the weighted (by priors Pi) likelihoods
      H(:,k)=log(const*Pi(k))+log_dM-(0.5*rsum(XM.*Xk));
      % we have estimated H(:,k) by its log!!
      logversions(k) = 1;
      disp(['Warning: Logversion for ' ks]);
    else
      %Xk is data - mean_k
      Xk=(data-ones(N,1)*Mu(:,k)'); 
      XM=Xk*MM;
      % H are the weighted (by priors Pi) likelihoods
      H(:,k)=const*Pi(k)*dM*exp(-0.5*rsum(XM.*Xk)); 	
    end
    
     components{k}.LocalFactors= XM*Lambda_k;
     
     Factors(:,factorCount+1:factorCount + numFactors(k))=components{k}.LocalFactors;
    
     factorCount = factorCount + numFactors(k);
  end
  
  Hsum=rsum(H);
  Hzero=(Hsum==0);  %the entries of H with zero sum
  Nz=sum(Hzero);  %the number of such entries
  H(Hzero,:)=eps*ones(Nz,numMeans)/numMeans;  %replace those with eps, equal probability
  Hsum(Hzero)=eps*ones(Nz,1);
  
  H=rdiv(H,Hsum);
  
 for k=1:numMeans
    components{k}.H=H(:,k);% responsibility : E[k|x];
 end

  