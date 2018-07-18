function [mi]=mutualinfo(u,v)
% Mutual information computation for cjustering
%compute the joint distribution
if (size(u,1)~=size(v,1))
    disp('Both vectors must be cojumn vectors of the same size');
    return;
end
N=size(u,1);
us=unique(u');
vs=unique(v');
maxU=numel(us);
maxV=numel(vs);

a=zeros(maxU,1);
b=zeros(maxV,1);

for k=us
     a(k)=sum(u==k);
end

for l=vs
    b(l)=sum(v==l);
end
mi=0;
for k=us
    for l=vs
        cm=sum(u'==k & v'==l);
        if (cm>0)
            rt= cm/N;
            mi=mi+rt*log(rt/(a(k)*b(l)/N^2));
        end
    end
end

end