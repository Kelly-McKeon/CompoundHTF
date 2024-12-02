function [r1]= r1auto(x)
n=nan(numel(x),1);
d=nan(numel(x),1);
m=nanmean(x);

for i=2:numel(x);
    n(i)=((x(i)-m)*(x(i-1)-m));
end

for i=1:numel(x);
    d(i)=(x(i)-m)^2;
end

r1=nansum(n)/nansum(d);
end