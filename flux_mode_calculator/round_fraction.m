function [A,l]=round_fraction(A_in,tol)
% Round floating point to nearest integer fraction with tolerance 'tol'.
% 'l' is the array with least common multiples of the denominator of the 
% rounded fractions, such that A*diag(l) is integer. 
%
%  Copyright (c) 2015, Jan Bert van Klinken
%

if nargin==1
    tol=1e-12;
end;

A=floor(A_in);
A_in=A_in-A;
if isempty(A_in)
    l=ones(1,0);
    return
end;

if nargout==1
    if ~issparse(A_in)
        [x,y]=calculate_fraction(A_in,tol);
        A=A+(x./y);
    else
        for i=find(full(sum(A_in~=0,1))~=0)
            a=full(A_in(:,i));
            j=a~=0;
            [x,y]=calculate_fraction(a(j,:),tol);
            A(j,i)=A(j,i)+sparse(x./y);
        end;
    end;
else
    if ~issparse(A_in)
        [x,y]=calculate_fraction(A_in,tol);
        A=A+(x./y);
        l=lcm_columns(y);
    else
        l=ones(1,size(A,2));
        for i=find(full(sum(A_in~=0,1))~=0)
            a=full(A_in(:,i));
            j=a~=0;
            [x,y]=calculate_fraction(a(j,:),tol);
            A(j,i)=A(j,i)+sparse(x./y);
            l(i)=lcm_columns(y(:));
        end;
    end;
    l((max(abs(A),[],1).*l)>=1e8)=Inf;
end;




