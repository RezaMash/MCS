function [X_out,x_left,x_right]=condition_matrix(X,varargin)
% Find coefficient vectors 'x_left' and 'x_right' such that the difference
% between the largest and smallest element in X_out=X.*(x_left*x_right) is 
% minimised. 
% condition_matrix(X,'Round2power','on') will round 'x_left' and 'x_right' 
% to the nearest powers of 2, such that the mantissa of the floating 
% points in 'X' will remain invariant.
%
%  Copyright (c) 2015, Jan Bert van Klinken
%

param=function_arguments(struct('Round2power','on'),varargin);


nz_left=full(sum(X~=0,2));
nz_right=full(sum(X~=0,1));

x_left=ones(size(X,1),1);
x_right=ones(1,size(X,2));
X1=abs(full(X));
X1(X1==0)=NaN;

i=1;
x=X1(:);
x(x==0)=[];
b=max(x)/min(x);
b1=100*b;
while (i<1000)&((b1-b)/b1>1e-6)
    x=exp(-mean(log([max(X1,[],1);min(X1,[],1)]),1));
    %x=(full(sum(X1,1))./nz_right).^a;
    x(nz_right==0)=1;
    x_right=x_right.*x;
    if issparse(X1)
        X1=X1*sparse(1:length(x),1:length(x),x);
    else
        X1=X1.*x(ones(size(X1,1),1),:);
    end;

    x=exp(-mean(log([max(X1,[],2),min(X1,[],2)]),2));
    %x=(full(sum(X1,2))./nz_left).^a;
    x(nz_left==0)=1;
    x_left=x_left.*x;
    if issparse(X1)
        X1=sparse(1:length(x),1:length(x),x)*X1;
    else
        X1=X1.*x(:,ones(1,size(X1,2)));
    end;
    i=i+1;
    b1=b;
    b=max(max(X1))/min(min(X1));
end;

if strcmp(param.Round2power,'on')
    x_left=2.^round(log(x_left)/log(2));
    x_right=2.^round(log(x_right)/log(2));
end;

if issparse(X)
    X_out=sparse(1:length(x_left),1:length(x_left),x_left)*X*sparse(1:length(x_right),1:length(x_right),x_right);
else
    X_out=X.*(x_left*x_right);
end;
