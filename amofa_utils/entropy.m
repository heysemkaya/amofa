function [entropy_val]=entropy(u)
% entropy computation for clustering
N=size(u,1);

entropy_val=0;
us=unique(u');
maxU=numel(us);
a=zeros(maxU,1);
for i=us
     a(i)=sum(u==i);
     rt=a(i)/N;
     entropy_val=entropy_val+rt*log(rt);
end

entropy_val=-entropy_val;

end