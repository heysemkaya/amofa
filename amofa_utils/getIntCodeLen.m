function [len]= getIntCodeLen(k)
% returns the integer code length in a recursive manner. 
% For details see Rissanen's "universal prior for integers" 
%c=2.865064
if (numel(k)>1)
    disp('call by one K at a time');
    return;
end

if k==1 
    len=1;
elseif k>1
    lk=ceil(log2(k));
    len=lk+getIntCodeLen(lk);
else
    len=0;
end

end