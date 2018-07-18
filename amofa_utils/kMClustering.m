function [M, v]= kMClustering(D,k,ctr)
% D: Dataset, k: number of components, 
% [optional] ctr: initial means
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
maxIter=3;
[m, n]=size(D);
mds = randperm(m); 
% randperm to select initial centers from among objects     

if nargin<3
    for i=1:k            
       ctr(i,:)=D(mds(i),:); % set random objects as initial centers then      	
    end
end
assgmnt=zeros(m,1);       

for iter=1:maxIter    
    d=dist(D,ctr');  %calc distance vector to current centers 
    [C,I]=min(d,[],2); % returns the indices and dist. of closest centers    
    if I==assgmnt         
        break;  % convergence case            
    else
        assgmnt=I; % set assignment as I        
    end
    for i=1:k            
           ctr(i,:)=mean(D(I==i,:),1); 
           %get mean of objects assigned to ith cluster
    end
end
M=ctr; % Matrix of k-means (kxn)
v=assgmnt; % assingments of objects to clusters 
end
