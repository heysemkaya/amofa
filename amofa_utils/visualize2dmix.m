function visualize2dmix(mixture,data,labels)
%visualize the mixture in 2D
colors = ['rgbkcm'];
colsize=size(colors,2);
numSamples=size(data,1);
%use the first two dimensions of the dataset
figure(1)
hold on
if nargin<3 
     for i = 1:numSamples
      plot(data(i,1),data(i,2),'.');
    end   
else 
    for i = 1:numSamples
      coli=mod(labels(i)-1,colsize)+1;  
      plot(data(i,1),data(i,2),[colors(coli),'.']);
    end
end

for i=1:numel(mixture)
  componentCovariance = mixture{i}.Lambda*mixture{i}.Lambda' + diag(mixture{i}.Psi);
  componentMean = mixture{i}.Mu;
  %use the first two dimensions 
  partialCov = componentCovariance(1:2,1:2)*4; %draw at 2 st. dev. line
  partialMean = componentMean(1:2);
  r = mvnplot2(partialMean,partialCov);
  plot(r(:,1),r(:,2),'-k');
end

end