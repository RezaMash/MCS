function [x_opt,f_opt,exitflag,output]=LP_solve(c,A,b,A_eq,b_eq,varargin)
% Solve Linear Programming problem:  
%
%     min c'*x subject to A*x <= b and A_eq*x = b_eq
%

param=function_arguments(struct('Solver','glpk','Tolerance',1e-10,'Integer','on','Bound',Inf,'ConditionCoeffs',{{}}),varargin);


if strcmp(lower(param.Integer),'on')
    param.Integer=(sum(sum(c~=round(c)))==0)&(sum(sum(A~=round(A)))==0)&(sum(sum(b~=round(b)))==0)&(sum(sum(A_eq~=round(A_eq)))==0)&(sum(sum(b_eq~=round(b_eq)))==0);
elseif strcmp(lower(param.Integer),'off')
    param.Integer=1==0;
else
    error(['Unknown option for ''Integer'': ',param.Integer])
end;

c=c(:);
m=[size(A,1),size(A_eq,1)];
n=length(c);
if isempty(A)&isempty(b)
    A=zeros(0,n);
    b=zeros(0,1);
end;
if isempty(A_eq)&isempty(b_eq)
    A_eq=zeros(0,n);
    b_eq=zeros(0,1);
end;
if (n~=size(A,2))|(n~=size(A_eq,2))
    error('Length of c should be equal to column size of A and A_eq . . .')
end;
if (m(1)~=size(A,1))
    error('Length of b should be equal to row size of A . . .')
end;
if (m(2)~=size(A_eq,1))
    error('Length of b_eq should be equal to row size of A_eq . . .')
end;
if n==0
    c=0;
    A_eq=zeros(size(A_eq,1),1);
    A=zeros(size(A,1),1);
    n=1;
end;    

if strcmp(lower(param.Solver),'glpk')
    % based on the MEX file of the OptiToolbox
    if exist('glpk')~=3
        disp(' ')
        disp('Error: Cannot find the MEX-function ''glpk'' . . .')
        disp(' ')
        disp('In order to install glpk: ')
        disp(' 1) get the glpkcc.cpp file from OptiToolbox and copy it to the current directory')
        disp(' 2) install glpk-448 and run Build_GLPK_with_VC9.bat')
        disp(' 3) give glpk_dir the full path name of the directory in which glpk-4.48 is stored')
        disp(' 4) eval([''mex -v -largeArrayDims glpkcc.cpp -I'',glpk_dir,''\src -L'',glpk_dir,''\w64 LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:libcmt.lib" -lglpk -lut -output glpk''])')
        disp(' ')
        error(' ')
    end;
    clear glpk;
    opts=struct;
    opts.msglev=1;
    a_in=full([A;A_eq]);
    b_in=full([b;b_eq]);
    c_in=full(c);
    if param.Integer
        % NB glpk_exact has a small memory leak and should not be run too often
        opts.lpsolver=3;
        opts.itlim=1000;   % prevent GLPK_exact from entering an infinite loop
        [x_opt1,f_opt1,exitflag1,output1]=glpk(c_in,a_in,b_in,repmat(-param.Bound,n,1),repmat(param.Bound,n,1),[repmat('U',m(1),1);repmat('S',m(2),1)],repmat('C',n,1),1,opts);
        if ~isinf(param.Bound)
            if sum(abs(abs(x_opt1)/param.Bound-1)<1e-8)>0
                exitflag1=6;
            end;
        end;
    else
        opts.lpsolver=1;
        [x_opt1,f_opt1,exitflag1,output1]=glpk(c_in,a_in,b_in,repmat(-Inf,n,1),repmat(Inf,n,1),[repmat('U',m(1),1);repmat('S',m(2),1)],repmat('C',n,1),1,opts);
    end;
    x_opt=x_opt1;
    f_opt=f_opt1;
    if exitflag1==5
        exitflag=1; % solution converged
    elseif (exitflag1==107)|(exitflag1==108)
        exitflag=0; % number of iterations exceeded
    elseif exitflag1==6
        exitflag=-2; % search space unbounded
    else
        exitflag=-1; % solution infeasible/error
    end;
    output=output1;
    
elseif strcmp(lower(param.Solver),'cplex')
    if ~exist('cplexlp')
        disp(' ')
        disp('Error: Cannot find the MEX-function ''cplexlp'' . . .')
        disp(' ')
        error(' ')
    end;
    [x_opt,f_opt,exitflag1,output]=cplexlp(c,A,b,A_eq,b_eq,repmat(-param.Bound,n,1),repmat(param.Bound,n,1),[]);
    
    exitflag=exitflag1>0;
    
else
    % add other solvers here
    error(['Unknown solver: ',param.Solver])
end;



if exitflag==1
    if ~isempty(param.ConditionCoeffs)
        s1a=1./param.ConditionCoeffs{2}(:);
        s1b=1./param.ConditionCoeffs{3}(:);
        s2=1./param.ConditionCoeffs{1}(:);
    else
        s1a=ones(size(A,1),1);
        s1b=ones(size(A_eq,1),1);
        s2=ones(size(A,2),1);
    end;
    i=(sum(A~=0,2)==1)&(b==0);
    A1=A(~i,:);
    b1=b(~i,:);
    s1a=s1a(~i,:);

    a=sum(abs(A1),1)./sum(A1~=0,1);
    a(a==0)=1;
    a(isnan(a))=1;
    x_opt(abs(x_opt.*a')./mean(abs(x_opt.*a'))<param.Tolerance)=0;
    i=abs(A1*x_opt-b1)./(abs(A1)*abs(x_opt)+abs(b1));
    output.A=[A1(isnan(i)|(i<param.Tolerance),:);A_eq];
    output.b=[b1(isnan(i)|(i<param.Tolerance),:);b_eq];
    s1=[s1a(isnan(i)|(i<param.Tolerance),:);s1b];
    i=x_opt==0;
    C=full(sparse(1:sum(i),find(i),1,sum(i),length(i)));
    output.A(:,i)=0;
    output.A=[output.A;C];
    output.b=[output.b;zeros(size(C,1),1)];
    s2(i)=1;
    s1=[s1;ones(sum(i),1)];

    output.A=output.A.*(s1*(s2'));
    output.b=output.b.*s1;
    

    if sum(((A*x_opt-b)./max(abs(A)*abs(x_opt)+abs(b),1e-100))>1e-9)>0
        % inequality constraints not satisfied
        exitflag=-1;
    end;
    if sum((abs(A_eq*x_opt-b_eq)./max(abs(A_eq)*abs(x_opt)+abs(b_eq),1e-100))>1e-9)>0
        % equality constraints not satisfied
        exitflag=-1;
    end;
end;

if exitflag~=1
    x_opt(:)=NaN;
    f_opt=NaN;
end;

