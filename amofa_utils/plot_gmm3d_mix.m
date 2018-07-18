function plot_gmm3d_mix(mixture,data,targets)
    d=size(data,2);
    if (d~=3)
       disp('Data must be numSamples x numDims with 3 dimensions');
       return;
    end
    usecolors=0;
    if nargin==3
        usecolors=1;
        colors = ['rgbcmyk'];
    end

    
    numSamples = size(data,1);
    numComponents = size(mixture,1);
    %use the first two dimensions of the iris dataset
    figure,
    hold on

    for i = 1:numSamples
        if (usecolors)
            plot3(data(i,1),data(i,2),data(i,3),[colors(mod(targets(i)-1,7)+1),'.']);% 
        else
            plot3(data(i,1),data(i,2),data(i,3),['b','.']);
        end
    end

    for i=1:numComponents
      componentCovariance = mixture{i}.Cov;
      componentMean = mixture{i}.Mu;
      %use the first two dimensions 
      partialCov = componentCovariance(1:3,1:3); 
      partialMean = componentMean(1:3);
      [hh] = plot_gaussian(partialCov, partialMean, 1, 10); %visualization function written by M.J. Beal
      %r = mvnplot2(partialMean,partialCov);
      %plot(r(:,1),r(:,2),'-k');
    end

   
    hold off
end