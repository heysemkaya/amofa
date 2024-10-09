# amofa
MATLAB scripts for AMoFA algorithm that is presented in https://arxiv.org/abs/1507.02801.
To run the code use 
>[mixture,training_history]=amofa(training_set);
where training_set is numSamples x numDimensions

For demo runs on several benchmark clustering and classification tasks see: 
demo_amofa.m

-----------------------------------------------------------------------------------------------------
This software is provided by the authors of the following paper for optimal transparency and comparison.
It is freely available for academic use given that the users should cite the relevant paper:

Heysem Kaya and Albert Ali Salah, Adaptive Mixtures of Factor Analyzers, submitted to Pattern Recognition

AMoFA is an automatic model selection algorithm for mixtures of factor analyzers. It finds the number of components and
their local latent dimensionality simultaneously with parameter estimation.
The algorithm is incremental-decremental. The incremental stage utilizes the fast component selection 
methods of Incremental Mixtures of Factor Analyzers (IMoFA) and improves its component splitting module. The proposed 
algorithm implemented here also contains IMoFA for comparison. 
The public code of IMoFA is also available from http://www.cmpe.boun.edu.tr/~salah/.

Unlike IMoFA, AMoFA does not resort to a validation set to tune the complexity. The proposed algorithm uses a
Bayesian-related MML criterion similar to the one proposed by Figueiredo and Jain in their paper titled 
"Unsupervised Learning of Finite Mixture Models" (The MATLAB code can be obtained from http://www.lx.it.pt/~mtf/).
The corresponding mixture code is given in folder ulfmm (with a slight modification, to avoid exhaustive plotting of models).
We generalize the MML criterion used there to reflect the local parameter cost as it should be balanced by the local support.
Components having less support than half of their parameter cost are automatically annihilated during the M-step of EM algorithm. 
Apart from this, the algorithm has a downsizing step that starts when the incremental stage is saturated. The downsizing stage annihilates
the weakest remaining component till the last component is left and runs EM iterations in between. 



For any questions/feedback and comments the corresponding author's email address: 
kaya.heysem@gmail.com

 
 Copyleft (2014): Heysem Kaya and Albert Ali Salah

 This software is distributed under the terms
 of the GNU General Public License Version 3
 
 Permission to use, copy, and distribute this software for
 any purpose without fee is hereby granted, provided that this entire
 notice is included in all copies of any software which is or includes
 a copy or modification of this software and in all copies of the
 supporting documentation for such software.
 This software is being provided "as is", without any express or
 implied warranty.  In particular, the authors do not make any
 representation or warranty of any kind concerning the merchantability
 of this software or its fitness for any particular purpose.








