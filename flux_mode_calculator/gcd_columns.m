function g=gcd_columns(X)
% calculates the column-wise greatest common divisors of matrix X
%
%  Copyright (c) 2015, Jan Bert van Klinken
%

g=zeros(1,size(X,2))+NaN;
i_integer=sum(X~=round(X),1)==0;
X=X(:,i_integer);

if ~isempty(X)
    if size(X,1)>1
        A=sort(abs(round(full(X))),1);
        A=[diff([zeros(1,size(A,2));A],[],1)>0].*A;
        a1=A(1,:);
        for i=2:size(A,1)
            b1=A(i,:);
            i1=b1>0;
            while sum(sum(i1))>0
                t=b1;
                b1(i1)=a1(i1)-floor(a1(i1)./b1(i1)).*b1(i1);
                a1(i1)=t(i1);
                i1=b1>0;
            end;
        end;
        c=a1+[a1==0];
    else
        c=abs(X)+[X==0];
    end;
    g(i_integer)=c;
end;
