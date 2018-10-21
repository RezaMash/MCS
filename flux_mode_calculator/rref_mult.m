function [X,Y]=rref_mult(A,B,C,D,tol)

if nargin<5
    tol=1e-10;
end;
if nargin<3
    C=[];
end;
if nargin<4
    D=[];
end;

if isempty(C)
    C=zeros(0,size(A,2));
end;
if isempty(D)
    D=zeros(size(C,1),size(B,2));
end;

Z=[A,B,zeros(size(B,1),1);C,D,ones(size(D,1),1)];

if sum(sum(Z~=round(Z)))==0
    [Z1,p,i]=rref_fast(Z,'ExcludeRows',[1:size(C,1)]+size(A,1),'ColumnOrder',1:size(A,2),'Integer','on','Tolerance',tol);
else
    [Z1,p,i]=rref_fast(Z,'ExcludeRows',[1:size(C,1)]+size(A,1),'ColumnOrder',1:size(A,2),'Tolerance',tol);
end;
if sum(sum(Z1~=round(Z1)))==0
    i(i<=size(A,1))=1:size(A,1);
    [a,i]=sort(i);
    Z1=Z1(i,:);
    X={Z1(1:size(A,2),size(A,2)+[1:size(B,2)]),diag(Z1(1:size(A,2),1:size(A,2)))};
    Y={-Z1(size(A,1)+[1:size(C,1)],size(A,2)+[1:size(B,2)]),Z1(size(A,1)+[1:size(C,1)],(size(A,2)+size(B,2)+1))};
    i=sign(Y{2});
    Y={Y{1}.*i,abs(Y{2})};
else
    i(i<=size(A,1))=1:size(A,1);
    [a,i]=sort(i);
    Z1=Z1(i,:);
    X=Z1(1:size(A,1),size(A,2)+[1:size(B,2)]);
    Y=-Z1(size(A,1)+[1:size(C,1)],size(A,2)+[1:size(B,2)])./Z1(size(A,1)+[1:size(C,1)],(size(A,2)+size(B,2)+1)*ones(1,size(B,2)));
    X(X~=0)=A(:,X~=0)\B;
    Y(Y~=0)=C(Y~=0,X~=0)*X(X~=0,:)+D(Y~=0,:);
end

if sum(p~=[ones(1,size(A,2)),zeros(1,size(B,2)+1)])>0
    disp('Warning: A not full column rank . . .')
end;

