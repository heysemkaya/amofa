function [Lambda,Mu,Psi,Pi,numFactors] = unpack(mixture);
%[Lambda,Mu,Psi,Pi,numFactors] = unpack(mixture);
% unpacks the mixture parameters into separate matrices
Lambda = [];
Mu = [];
Psi = [];
Pi = [];
numFactors = [];
for i = 1:size(mixture,2)
  Lambda = [Lambda mixture{i}.Lambda];
  Mu = [Mu mixture{i}.Mu];
  Psi = [Psi mixture{i}.Psi];
  Pi = [Pi mixture{i}.Pi];
  numFactors = [numFactors mixture{i}.numFactors];
end