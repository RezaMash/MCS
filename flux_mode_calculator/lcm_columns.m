function c=lcm_columns(X)
% calculates the column-wise least common multiples of matrix X
%
%  Copyright (c) 2015, Jan Bert van Klinken
%

A=sort(abs(round(full(X))),1);
A=[diff([zeros(1,size(A,2));A],[],1)>0].*A;
A(A==0)=1;
pA=prod(A,1);
c=pA./gcd_columns([ones(size(A,1),1)*pA]./A);