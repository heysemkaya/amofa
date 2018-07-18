function plot_mixture(mixture,data,targets)
    usecolors=0;
    if nargin==3
        usecolors=1;
        colors = ['rgcmykb'];
    end

    
    numSamples = size(data,1);
    numComponents = size(mixture,1);

    hold on
    for i = 1:numSamples
        if (usecolors)
            plot(data(i,1),data(i,2),[colors(mod(targets(i)-1,7)+1),'.']);% 
        else
            plot(data(i,1),data(i,2),['c','.']);
        end
    end

    for i=1:numComponents
      componentCovariance = mixture{i}.Lambda*mixture{i}.Lambda' + diag(mixture{i}.Psi);
      componentMean = mixture{i}.Mu;
      %use the first two dimensions 
      partialCov = componentCovariance(1:2,1:2)*4; %draw at 2 st. dev. line
      partialMean = componentMean(1:2);
      r = mvnplot2(partialMean,partialCov);
      plot(r(:,1),r(:,2),'-k');
    end
    hold off
end