function plot_mixture3d(mixture,data,targets)
    usecolors=0;
    if nargin==3
        usecolors=1;
        colors = ['rgbcmyk'];
    end

    
    numSamples = size(data,1);
    numComponents = size(mixture,1);
    %use the first two dimensions of the iris dataset
    %figure,
    hold on

    %for i = 1:numSamples
        if (usecolors)
            %plot3(data(i,1),data(i,2),data(i,3),[colors(mod(targets(i)-1,7)+1),'.']);% 
            scatter3(data(:,1),data(:,2),data(:,3),[colors(mod(targets-1,7)+1),'.']);
        else
            scatter3(data(:,1),data(:,2),data(:,3),['b','.']);
        end
    %end

    for i=1:numComponents
      componentCovariance = mixture{i}.Lambda*mixture{i}.Lambda' + diag(mixture{i}.Psi);
      componentMean = mixture{i}.Mu;
      %use the first two dimensions 
      partialCov = componentCovariance(1:3,1:3)*4; 
      partialMean = componentMean(1:3);
      [hh] = plot_gaussian(partialCov, partialMean, 1, 10); %visualization function written by M.J. Beal

    end

   
    hold off
end