function visualize2d(bestk,bestmu,bestcov,data,labels)
%visualize the mixture in 2D
colors = ['rgbcmyk'];
colsize=size(colors,2);
numSamples=size(data,1);
%use the first two dimensions of data
figure(1)
hold on
if nargin<5 
     for i = 1:numSamples
      plot(data(i,1),data(i,2),'.');
    end   
else 
    for i = 1:numSamples
      coli=mod(labels(i)-1,colsize)+1;  
      plot(data(i,1),data(i,2),[colors(coli),'.']);
    end
end
for i=1:bestk
  componentCovariance = bestcov(:,:,i);
  componentMean = bestmu(:,i);
  %use the first two dimensions 
  partialCov = componentCovariance(1:2,1:2)*4; %draw at 2 st. dev. line
  partialMean = componentMean(1:2);
  r = mvnplot2(partialMean,partialCov);
  plot(r(:,1),r(:,2),'-k');
end

end