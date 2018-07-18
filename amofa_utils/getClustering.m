function  [clustering,responsibility,maxlik,params]= getClustering(data,mixture)
% returns the hard labeling and the soft responsibilities of given data for
% given mixture 
% clustering: hard labeling of points with respect to cluster indices Nx1
% responsibility: normalized component responsibilities NxK
% maxlik: maximum normalized log likelihood for each 


      params=0;
      [Lambda,Mu,Psi,Pi,numFactors] = unpack(mixture);
      % data,lambda,psi,mu,priors,numFactors,regTerm
      loglikeVal = loglike4(data,Lambda,Psi,Mu,Pi,numFactors);
      sumLoglik=logsumexp(loglikeVal,2);
      loglikeVal=loglikeVal-repmat(sumLoglik,1,size(loglikeVal,2));
      params=params+ numel(Lambda)+numel(Psi)+numel(Mu)+numel(Pi);
      responsibility=exp(loglikeVal);
      %if size(loglikeVal,2)>1
      %more than one component
      [maxlik, clustering] = max(loglikeVal,[],2);

end