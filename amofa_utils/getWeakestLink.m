% -----------------------------------------------------------------------
% Copyleft (2014): Heysem Kaya and Albert Ali Salah
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

function [weakestLink] =getWeakestLink(Pi,numSamples,numDims,numFactors,MinSoftCount)
% [weakestLink] =getWeakestLink(Pi,numSamples,numDims,numFactors,MinSoftCount)
% Selects the weakest component whose support is less than
% parametric/theoretic threshold

    weakestLink=0;
    if (MinSoftCount>0)
      [minVal,minIndex] = min(Pi*numSamples);
      if (minVal<MinSoftCount)
       weakestLink=minIndex;
      end
    elseif (MinSoftCount==-1)
        %compute weakness using theoretic value of C_j
        %for all components individually
        ComponentParams=numDims*(numFactors+2)+2; %factors + means + psi + 1 ( pi(k))
        [minVal,minIndex]=min(Pi'*numSamples-(ComponentParams/2));
        if (minVal<=0)
            weakestLink=minIndex;
        end
    end
end