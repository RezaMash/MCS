function c=gcd_fast(a,b)
% vectorised computation of the greatest common divisor of (matrix) a and b

b1=abs(round(b));
a1=abs(round(a));
i2=(abs(b-b1)>1e-12)|(abs(a-a1)>1e-12);
i1=b1>0;
while sum(sum(i1))>0
    t=b1;
    b1(i1)=a1(i1)-floor(a1(i1)./b1(i1)).*b1(i1);
    a1(i1)=t(i1);
    i1=b1>0;
end;
c=a1+[a1==0];
c(i2)=NaN;
