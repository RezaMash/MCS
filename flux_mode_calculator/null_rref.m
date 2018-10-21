function N=null_rref(A,pivot)
% calculate the kernel matrix 'N' of matrix 'A', where 'A' is in reduced 
% row echelon form and 'pivot' is a logical array with the pivots of 'A' 
%
%  Copyright (c) 2015, Jan Bert van Klinken
%

if sum(pivot)==0
    if issparse(A)
        N=sparse(1:length(pivot),1:length(pivot),1);
    else
        N=eye(length(pivot));
    end;
    return
end;

pivot=pivot~=0;

if isnumeric(A)
    i=sum(A(:,pivot)~=0,2);
    if sum(i>1)>0
        error('Input matrix is not in reduced row echelon form . . .') 
    end;
    A=A(i==1,:);
    D=A(:,pivot);
    n=size(D,1);
    if [sum(sum(D~=0,1)~=0)~=n]|[sum(sum(D~=0,2)~=0)~=n]
        error('Input matrix is not in reduced row echelon form . . .') 
    end;
    p=[1:n]*[D~=0];
    d=diag(D(p,:));
    if sum(sum(A~=round(A)))==0
        l=lcm_columns(d(:));
        if (isnan(l))|(l>=1e6)
            l=1;
        end;
    else
        l=1;
    end;
    D_inv=sparse(1:n,p,l./d);
    if issparse(A)
        N=sparse(size(A,2),sum(~pivot));
        N(pivot,:)=-D_inv*A(1:sum(pivot),~pivot);
        N(~pivot,:)=sparse(1:sum(~pivot),1:sum(~pivot),l);
    else
        N=zeros(size(A,2),sum(~pivot));
        N(pivot,:)=-D_inv*A(1:sum(pivot),~pivot);
        N(~pivot,:)=l*eye(sum(~pivot));
    end;
    if sum(sum(A~=round(A)))==0
        x=sum(N~=round(N),1)>0;
        
        % divide by gcd for integers
        g=gcd_columns(full(N(:,~x)));
        N(:,~x)=round(N(:,~x)./g(ones(size(N,1),1),:));

        % perform fraction rounding for floats
        x=find(x);
        [a,l]=round_fraction(N(:,x),1e-12);
        x=x(l<1e6);
        l=l(:,l<1e6);
        N(:,x)=round(N(:,x).*l(ones(size(N,1),1),:));        
    end
else
    error('A should be numeric . . .')
end;



