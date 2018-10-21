function [A,pivot,p_row,row_indep]=rref_fast(A_in,varargin)
% Compute reduced row echelon form of A_in. If A_in is non-integer, left 
% division is used to determine linearly dependent rows and columns, to
% detect true zeros, and to increase numerical precision.
%
%  Copyright (c) 2015, Jan Bert van Klinken
%

param=function_arguments(struct('ColumnOrder',[],'ExcludeRows',[],'Tolerance',1e-10,'PivotSelection','partial','NonNegativityConstraints',[],'Pivots',zeros(1,size(A_in,2)),'ReferenceMatrix',A_in,'Integer','off'),varargin);



if ~islogical(param.ExcludeRows)
    param.ExcludeRows=full(sparse(1,param.ExcludeRows,1,1,size(A_in,1)))~=0;
end;

if isempty(A_in(~param.ExcludeRows,:))
    A=A_in;
    pivot=zeros(1,size(A_in,2))==1;
    p_row=[1:size(A,1)]';
    row_indep=[];
    
    return
end;

param.PivotSelection=lower(param.PivotSelection);
if sum(strcmp(param.PivotSelection,{'partial','full','partial compression','full compression'}))==0
    error(['Unknown type of pivot selection: "',param.PivotSelection,'" . . .'])
end;
if (sum(strcmp(param.PivotSelection,{'partial','full'}))==0)&(length(param.NonNegativityConstraints)~=size(A_in,2))
    error('NonNegativityConstraints argument should be given in compression mode . . .')
end;

if ~isempty(param.ColumnOrder)
    param.ColumnOrder=param.ColumnOrder(:)';
    param.ColumnOrder=[param.ColumnOrder(:,param.Pivots(param.ColumnOrder)==1),param.ColumnOrder(:,param.Pivots(param.ColumnOrder)==0)];
    i_column=[param.ColumnOrder,find(full(sparse(1,param.ColumnOrder,1,1,size(A_in,2)))==0)];
    k1=length(param.ColumnOrder);
    A_in=A_in(:,i_column);
    param.ReferenceMatrix=param.ReferenceMatrix(:,i_column);
    param.Pivots=param.Pivots(:,i_column);
    if ~isempty(param.NonNegativityConstraints)
        param.NonNegativityConstraints=param.NonNegativityConstraints(i_column);
    end;
else
    k1=size(A_in,2);
end;

if ~isnumeric(A_in)
    error('Input matrix should be numeric . . .')
end;

if strcmp(param.Integer,'on')
    param.all_integer=sum(sum(A_in~=round(A_in)))==0;
elseif strcmp(param.Integer,'off')
    param.all_integer=1==0;
else
    error(['Unknown option for ''Integer'': ',param.Integer,' . . .'])
end;

% order rows such that excluded rows remain last
p_row=[find(~param.ExcludeRows(:));find(param.ExcludeRows(:))];
A=A_in(p_row,:);
A_in=param.ReferenceMatrix;

% ensure first rows are pivots
pivot=param.Pivots*2;
a=double(A(1:sum(pivot==2),pivot==2)~=0);
if ~isempty(a)
    a=a*a'-eye(size(a,1));
    if sum(sum(a~=0))>0
        error('First rows are not pivot rows . . .')
    end;
end;


% calculate rank
[A_in,a_left,a_right]=condition_matrix(A_in,'Round2power','on');
x=sort(svd(full(A_in(~param.ExcludeRows,:))));
if sum(x>1e-30)>0
    x(x==0)=min(1e-16,min(x(x~=0)));
    k=exp(diff(log(x)));
    k(isnan(k)|isinf(k))=1;
    k=[0;find(k>1e9)]+1;
    k=length(x)-k(end)+1;
else
    k=0;
end;

[n,m]=size(A);

rev=1-param.NonNegativityConstraints(:)';


j=sum(pivot==2)+1;
i=-1;
nz=-1;

while (sum(pivot~=0)<k1)&(sum(pivot==2)<k)&(strcmp(param.PivotSelection,'partial')|(~isempty(nz)))
    if strcmp(param.PivotSelection,'partial')
        i=find(pivot==0);
        i=i(1);
        nz=find((A(:,i)'~=0)&([1:n]>=j)&(~param.ExcludeRows(:,p_row)));
        
    elseif strcmp(param.PivotSelection,'full')
        I=zeros(sum(sum(A~=0)),3);
        [I(:,1),I(:,2),I(:,3)]=find(A);
        I(pivot(I(:,2))~=0,:)=[];
        I(I(:,1)<j,:)=[];
        I(I(:,2)>k1,:)=[];
        I(param.ExcludeRows(:,p_row(I(:,1))),:)=[];
        
        [a,i]=sort(-abs(I(:,3)));
        I=I(i,:);
        
        if ~isempty(I)
            nz=I(1,1);
            i=I(1,2);
        else
            i=0;
            nz=[];
        end;

    else
        % partial or full compression
        N=zeros(size(A,1),4);
        I=zeros(size(A,1),4);
        nz=[1:size(A,1)]';
        
        B=(A<0).*(1-rev(ones(size(A,1),1),:));
        N(:,1)=sum(B,2);
        [a,I(:,1)]=max(abs(A).*B,[],2);
        B=(A>0).*(1-rev(ones(size(A,1),1),:));
        N(:,2)=sum(B,2);
        [a,I(:,2)]=max(abs(A).*B,[],2);
        B=(A~=0).*(rev(ones(size(A,1),1),:));
        N(:,3)=sum(B,2);
        [a,I(:,3)]=max(abs(A).*B,[],2);

        if strcmp(param.PivotSelection,'full compression')
            t=1*((N(:,1)==1)&(N(:,2)~=0)&(N(:,3)==0)) + 2*((N(:,1)~=0)&(N(:,2)==1)&(N(:,3)==0)) + 4*((N(:,1)==0)&(N(:,2)~=0)&(N(:,3)==1)) + 4*((N(:,1)~=0)&(N(:,2)==0)&(N(:,3)==1)) + 5*((N(:,1)==0)&(N(:,2)==0)&(N(:,3)==2));
        else
            % only include coupled reaction pairs
            t=3*((N(:,1)==1)&(N(:,2)==1)&(N(:,3)==0)) + 4*((N(:,1)==0)&(N(:,2)==1)&(N(:,3)==1)) + 4*((N(:,1)==1)&(N(:,2)==0)&(N(:,3)==1)) + 5*((N(:,1)==0)&(N(:,2)==0)&(N(:,3)==2));
        end;
        
        J=nz*0+1;
        J(t==1,:)=I(t==1,1);
        J(t==2,:)=I(t==2,2);
        J(t==3,:)=I(t==3,2);
        J(t>=4,:)=I(t>=4,3);

        i=(nz<j)|(param.ExcludeRows(:,p_row(nz))')|(t==0)|(J>k1)|(pivot(:,J)'~=0);
        I(i,:)=[];
        J(i,:)=[];
        N(i,:)=[];
        nz(i,:)=[];
        t(i,:)=[];
        
        i=0;
        ii=find((sum(A~=0,2)==1)&(~param.ExcludeRows(:,p_row)')&([1:size(A,1)]'>sum(pivot==2)));
        ii(sum(A(ii,1:k1),2)==0)=[];

        if ~isempty(ii)
            nz=ii(1);
            i=find(A(nz,:)~=0);
        
        elseif ~isempty(nz)
            [J,i]=sort(J);
            nz=nz(i,:);
            i=J(1);
            nz=nz(1);
        end;


    end;

    if i>0
        if max(abs([NaN;A(nz,i)]))<0.01
            ii=p_row([1:sum(pivot==2),nz(:)']);
            x=A_in(ii,pivot==2)\A_in(ii,i);
            if mean(abs(A_in(ii,pivot==2)*x-A_in(ii,i)))/mean(abs(A_in(ii,pivot==2))*abs(x)+abs(A_in(ii,i)))<param.Tolerance
                x=correct_zeros(x,A_in(ii,pivot==2),A_in(ii,i),param.Tolerance);
                x=x.*((a_right(:,pivot==2)')*(1./a_right(:,i)));
                A(1:sum(pivot==2),i)=A(1:sum(pivot==2),pivot==2)*x;
                A(nz,i)=0;
                nz=[];
            end;
        end;
    end;
    
    if ~isempty(nz)
        % swap rows with row with 
        if param.all_integer
            % smallest nonzero element
            [b,jj]=min(abs(A(nz,i)));
        else
            % largest nonzero element
            [b,jj]=max(abs(A(nz,i)));
        end;
        ii=nz(jj);
        if ii~=j
            A([ii,j],:)=A([j,ii],:);
            p_row([ii,j])=p_row([j,ii]);
        end;

        i1=find((A(:,i)'~=0)&([1:n]~=j));
        i2=find(pivot==0);
        aa=full(A(j,i));
        if aa<0
            aa=-aa;
            A(j,:)=-A(j,:);
        end;

        if  param.all_integer
            A1=A(i1,:)*aa-A(i1,i)*A(j,:);            
            if full(max(max(abs(A1))))>2^53
                param.all_integer=1==0;
                a=sum(A(1:sum(pivot==2),pivot==2),2);
                a=1./a(:);
                if issparse(A)
                    A(1:sum(pivot==2),:)=sparse(1:sum(pivot==2),1:sum(pivot==2),a)*A(1:sum(pivot==2),:);
                else
                    A(1:sum(pivot==2),:)=a(:,ones(1,size(A,2))).*A(1:sum(pivot==2),:);
                end;
                A(j,i2)=A(j,i2)/aa;
                A1=A(i1,i2)-A(i1,i)*A(j,i2);
            end;
        else
            A(j,i2)=A(j,i2)/aa;
            A1=A(i1,i2)-A(i1,i)*A(j,i2);
        end;

        if param.all_integer
            g=gcd_columns(abs(A1'));
            g=g(:);
            ii=g>1;
            if issparse(A1)
                A1(ii,:)=sparse(1:sum(ii),1:sum(ii),1./g(ii))*A1(ii,:);
            else
                A1(ii,:)=A1(ii,:)./g(ii,ones(1,size(A1,2)));
            end;
            A(i1,:)=round(A1);
        else
            A(i1,i)=0;
            A(i1,i2)=A1;
        end;


        pivot(1,i)=2;
        j=j+1;
    elseif (i>0)&(sum(A((~param.ExcludeRows(:,p_row))&([1:size(A,1)]>sum(pivot==2)),i)~=0)==0)
        pivot(1,i)=1;
    end;
end;

if sum(param.ExcludeRows(p_row(1:sum(pivot==2))))>0
    error('ExcludeRows parameter was not satisfied . . .')
end;

% postprocessing
if param.all_integer
    a=gcd_columns(A')';
    a(isnan(a))=1;
    a=1./a;
    a(1:sum(pivot==2))=a(1:sum(pivot==2),:).*sign(sum(A(1:sum(pivot==2),pivot==2),2));
    if issparse(A)
        A=sparse(1:length(a),1:length(a),a)*A;
    else
        A=a(:,ones(1,size(A,2))).*A;
    end;
    
elseif sum(pivot==2)>0
    % pivot rows i
    i=1:sum(pivot==2);
    % rows j linearly dependent on the pivot rows
    j=(sum(pivot==2)+1):size(A_in,1);
    j=j(:,((sum(abs((A_in(p_row(i),:)')*((A_in(p_row(i),:)')\(A_in(p_row(j),:)'))-A_in(p_row(j),:)'),1)./sum(abs(A_in(p_row(j),:)')))<param.Tolerance)|(sum(A_in(p_row(j),:)'~=0,1)==0));

    % complement rows l of i and j
    l=1:size(A_in,1);
    l([i,j])=[];
    
    % check rank
    if rank(full(A_in(p_row([i,j]),pivot==2)))<sum(pivot==2)
        error('Invalid reference matrix . . .')
    end;
    x=A_in(p_row([i,j]),pivot==2)\A_in(p_row([i,j]),pivot~=2);
    x=correct_zeros(x,A_in(p_row([i,j]),pivot==2),A_in(p_row([i,j]),pivot~=2),param.Tolerance);
    A(i,pivot~=2)=A(1:sum(pivot==2),pivot==2)*(x.*((a_right(:,pivot==2)')*(1./a_right(:,pivot~=2))));
    A(j,:)=0;
    
    if ~isempty(l)
        A1=sparse(p_row(l),1:length(l),1,size(A_in,1),length(l));
        if ~issparse(A_in)
            A1=full(A1);
        end;

        
        y=zeros(length(l),sum(pivot~=2));
        yy=[zeros(1,sum(pivot~=2));x];
        for ii=1:length(l)
            jj=full(sum(A1(:,[1:length(l)]~=ii),2))==0;
            yy(1,:)=A_in(p_row(l(ii)),pivot~=2)-(A_in(p_row(l(ii)),pivot==2)*x);
            yy1=correct_zeros(yy,[A1(jj,ii),A_in(jj,pivot==2)],A_in(jj,pivot~=2),param.Tolerance);
            y(ii,:)=(yy1(1,:)./a_right(:,pivot~=2))/a_left(p_row(l(ii)));
        end;

        A(l,pivot~=2)=y;
        
    end;
end;


        
pivot=pivot==2;

if nargout>=4
    row_indep=full(sum(A(:,pivot),2))>0;
    j=find((~row_indep')&(sum(A_in(p_row,:)'~=0,1)>0)&(sum(A'~=0,1)>0));
    if ~isempty(j)
        A1=condition_matrix(A_in);
        for i=j
            row_indep(i)=(sum(abs((A1(p_row(find(row_indep(:)')),:)')*((A1(p_row(find(row_indep(:)')),:)')\(A1(p_row(i),:)'))-A1(p_row(i),:)'),1)./sum(abs(A1(p_row(i),:)')))>param.Tolerance;
        end;
    end;
end;

% validity check of compression results
if sum(strcmp(param.PivotSelection,{'partial','full'}))==0
    p=pivot;
    p(param.Pivots~=0)=0==1;
    if (sum(full(sum(A((1:sum(p))+sum(param.Pivots),~p)~=0,2))>1)>0)&strcmp(param.PivotSelection,'partial compression')
        error('Compression is not partial . . .');
    end;
    x=A((1:sum(p))+sum(param.Pivots),p&(rev==0))'*A((1:sum(p))+sum(param.Pivots),~p);
    if full(sum(sum(x>0)))>0
        error('Compression violated the non-negativity constraints . . .');
    end;
end;


if ~isempty(param.ColumnOrder)
    [a,i]=sort(i_column);
    A=A(:,i);
    A_in=A_in(:,i);
    pivot=pivot(:,i);
end;




function x=correct_zeros(x_in,A,B,tol)

x=x_in;
Y=full(B);
for i=find(sum(x~=0,1)>0)
    [y,j]=sort(abs(x(:,i)),1,'descend');
    k=sum(y>0.01);
    j(y==0)=[];
    
    X=full(A(:,j));    
    y=X(:,1:k)\Y(:,i);
    b=1/mean(abs(X(:,1:k))*abs(y)+abs(Y(:,i)));
    while ((mean(abs(X(:,1:k)*y-Y(:,i)))*b)>tol)&(k<length(j))
        k=k+1;
        y=X(:,1:k)\Y(:,i);
    end;

    x(j((k+1):end),i)=0;
    x(j(1:k),i)=y;
end;




