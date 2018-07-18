function [res]=norm_info_dist(u,v)

res= 1- (mutualinfo(u,v)/max([entropy(u) entropy(v)]));

end