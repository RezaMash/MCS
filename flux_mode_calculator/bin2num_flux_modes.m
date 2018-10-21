function [fm,err]=bin2num_flux_modes(fm_bin,S)
% [fm,err] = BIN2NUM_FLUX_MODES(fm_bin,S)
%  This function calculates the coefficients of the binary/ternary flux
%  modes (fm_bin) such that S*fm=0 and fm_bin=sign(fm), with S the
%  stoichiometric matrix. Importantly, if fm_bin represent intermediate
%  flux modes for which not all compounds have been constrained to zero
%  steady state, only the constrained rows of S must be included. Indices
%  of the constrained compounds are given in the stats.id variable that is
%  returned by the CALCULATE_FLUX_MODES function. The indices of
%  inconsistent binary flux modes are returned in the cell array err, with
%  in the first column the index of the flux mode and in the second column
%  the error messsage.
%
%
%  Copyright (c) 2015, Jan Bert van Klinken
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%    
%

tol=1e-10;
if size(S,2)~=size(fm_bin,1)
    error('Number of rows in fm_bin must be equal to the number of columns in S')
end;
S_in=S;
[a,l]=round_fraction(S,1e-12);
l=find((l>1e6)|full(max(abs(S),[],1)>1e5));
i_frac=ones(1,size(fm_bin,2))==1;
for i=1:length(l)
    i_frac=i_frac&(fm_bin(l(i),:)==0);
end;

if size(fm_bin,2)>size(fm_bin,1)
    [s,i]=sort(sum(S~=0,1));
    if sum(sum(S~=round(S)))>0
        [S,s_left,s_right]=condition_matrix(S);
        s_right=s_right(:);
        [S1,S_pivot]=rref_fast(S,'ColumnOrder',i,'Tolerance',tol,'Integer','off');
    else
        [S1,S_pivot]=rref_fast(S,'ColumnOrder',i,'Tolerance',tol,'Integer','on');
        [S1,s_left,s_right]=condition_matrix(S1);
        s_right=s_right(:);
        S1(1:sum(S_pivot),:)=S1(1:sum(S_pivot),:)./repmat(sum(S1(1:sum(S_pivot),S_pivot),2),1,size(S1,2));
    end;
    S=S1(1:sum(S_pivot),S_pivot)'*S1(1:sum(S_pivot),:);
    clear S1
    
    S=full(S');
    S_pivot=S_pivot(:);
else
    [S,s_left,s_right]=condition_matrix(S);
    S=full(S');
    s_right=s_right';
    S_pivot=[];
end;

fm=zeros(size(fm_bin));
err=cell(0,2);
X=zeros(size(fm_bin,1),1);
for i=1:size(fm_bin,2)
    err1='';
    a=double(fm_bin(:,i));
    if sum(a==-128)==0
        z=a~=0;
        a1=a(z,:);
        if ~isempty(S_pivot)
            A=S((~S_pivot)&z,~z(S_pivot));
        else
            A=S(z,:);
        end;
        
        if size(A,1)>size(A,2)
            A=[A,zeros(size(A,1),size(A,1)-size(A,2))];
        end;
        [U,s]=svd(A,0);
        if (length(a1)>1)&(size(s,1)>1)
            s=diag(s);
        elseif size(s,1)==1
            s=s(1);
        else
            s=[];
        end;
        if length(s)>1
            s1=sort(s);
            s1=s1(ceil((length(s1)-1)/2)+1);
            s1=s1+(s1<1e-10);
            s=s/s1;
        end;
        NS=U(:,s<tol);

        if size(NS,2)==1
            if ~isempty(S_pivot)
                X(S_pivot(z),:)=-S((~S_pivot)&z,z(S_pivot))'*NS;
                X(~S_pivot(z),:)=NS;
                X(abs(X)<tol)=0;
                NS=X(1:sum(z),:);
            end;
            NS=NS.*s_right(z,:);
            NS=NS*sign(NS(1))*sign(a1(1));
            if sum(sign(NS)~=sign(a1))==0             
                if i_frac(i)
                    NS=NS/min(abs(NS));
                    [NS,l]=round_fraction(NS,1e-12);
                    if l<1e6
                        NS=round(NS*l);    
                    end;
                end;
                fm(z,i)=NS;
            else
                err1=['Nullspace has dimension 1 but the signs of the flux mode are incorrect . . .'];
            end;
        elseif size(NS,2)>1
            err1=['Nullspace has dimension ',num2str(size(NS,2)),' . . .'];
        elseif size(NS,2)==0
            err1=['No nullspace found: minimal eigenvalue is ',num2str(min(s)),' . . .'];
        end;    
    else
        err1='Inconsistent flux mode due to unpacking . . .';
    end;
    if ~isempty(err1)
        err=[err;{i,err1}];
        fm(:,i)=NaN;
    end;
end;


