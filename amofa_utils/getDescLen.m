function [descLength,modelLogLik, numparams]=getDescLen(model,modelLogLik,ModelSelCri)
%Gets the description lenght in terms of MML (ModelSelCri=1), MDL/BIC
%(ModelSelCri=2), or AIC given the MoFA model and its log likelihood
% -----------------------------------------------------------------------
% Copyleft (2014): Heysem Kaya 
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

Mu=model.Mu;
numSamples=model.numSamples;
Lambda=model.Lambda;
Pi=model.Pi;
Psi=model.Psi;
numFactors=model.numFactors;
[numDims,numMeans]=size(Mu);
numparams=zeros(numMeans,1);
numFactors_bits=zeros(numMeans,1);
logc=log2(2.865064);
%Rissanen J., Information and Complexity in Statistical Modeling, Springer, 2007, p.15
for k=1:numMeans
    numparams(k)=numDims*(numFactors(k)+2)+1;
    numFactors_bits(k)= getIntCodeLen(numFactors(k))+logc;
end
%We keep a matrix for independent Psi however !

nparsover2= numparams/2;

code_len_k= getIntCodeLen(numMeans)+logc;

if (ModelSelCri==1) 
     descLength = -modelLogLik*log2(exp(1))+ (sum(nparsover2.*log2(Pi))) + ...
         sum(nparsover2 + 0.5)*log2(numSamples/12)+ ...
          code_len_k+sum(numFactors_bits); 
elseif (ModelSelCri==2)
        descLength=-modelLogLik + sum(nparsover2)*log(numSamples);
elseif  (ModelSelCri==3) 
     descLength=-modelLogLik + sum(nparsover2);
end
                
end