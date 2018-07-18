function plot_gmm(data,cov,mu,targets)
    d=size(data,2);
    if (d~=2)
       disp('Data must be numSamples x numDims with 2 dimensions');
       return;
    end
    usecolors=0;
    if nargin==4
        usecolors=1;
        colors = ['rgcmykb'];
    end

    
    numSamples = size(data,1);
    numComponents = size(mu,2);
    %use the first two dimensions of the iris dataset
    figure,
    hold on

    for i = 1:numSamples
        if (usecolors)
            plot(data(i,1),data(i,2),[colors(mod(targets(i)-1,7)+1),'.']);% 
        else
            plot(data(i,1),data(i,2),['b','.']);
        end
    end

    for i=1:numComponents
      componentCovariance = cov(:,:,i)*4;% 2 -std line
      componentMean = mu(:,i);

      r = mvnplot2(componentMean ,componentCovariance);
      plot(r(:,1),r(:,2),'-k');

    end

   
    hold off
end