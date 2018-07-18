function  [clustering,responsibility,maxlik,params]= getClusteringGmm(data,Mu,covM,Pi)
% returns the hard labeling and the soft responsibilities of given data for
% given gaussian mixture 
% clustering: hard labeling of points with respect to cluster indices Nx1
% responsibility: normalized component responsibilities NxK
% maxlik: maximum normalized log likelihood for each 

      params=0;
      % data,lambda,psi,mu,priors,numFactors,regTerm
      loglikeVal = loglike_gmm(data,Mu,covM,Pi,1e-4);
      sumLoglik=logsumexp(loglikeVal,2);
      loglikeVal=loglikeVal-repmat(sumLoglik,1,size(loglikeVal,2));
      params=params+ numel(covM)+numel(Mu)+numel(Pi);
      responsibility=exp(loglikeVal);
      %if size(loglikeVal,2)>1
      %more than one component
      [maxlik, clustering] = max(loglikeVal,[],2);

end
