function [r]=logsumexp(L,dim)
% Numerically stable computation of log(sum(exp(L))
% Heysem Kaya - Usage proposed by ATC
    if (nargin<2)
        dim=1;
    end
    
    mx=max(L,[],dim);
   
    
    %d(dim)=1;
    M=bsxfun(@minus, L,mx);
    r=mx+log(sum(exp(M),dim));

end