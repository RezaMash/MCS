function [fm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(S,c,varargin)
% [fm_bin,S_unc,id_unc,T,stats] = CALCULATE_FLUX_MODES(S,c,'PropertyName',PropertyValue,...)
%  This function calculates the elementary flux modes, extreme pathways, or
%  minimal generating sets of the stoichiometric matrix S with constraints
%  c, where c(i) = 1 for reversible reactions, c(i) = 0 for irreversible
%  reactions, and c(i) = -1 for zero (blocked) reactions. As output the
%  binary steady state flux modes fm_bin are returned and the coefficients
%  (S_unc) and identifiers (id_unc) of the unconstrained reactions and/or
%  compounds. Unconstrained reactions/compounds include pre-defined ones
%  using the 'ReactionsUnconstrained' and 'CompoundsUnconstrained' option,
%  and the remaining reactions/compounds on which no constraint could be
%  imposed because the number of intermediate flux modes surpassed
%  'MaxIntermediates' or because not enough memory was available to start a
%  new iteration. The first row of id_unc indicates whether the identifier
%  refers to a compound (1) or reaction (2), and the second row gives the
%  compound/reaction index. If the compression option is used, T will
%  return the transformation matrix such that the compressed stoichiometric
%  matrix is equal to S*T. Importantly, the reaction indices in id_unc
%  refer to the columns in the compressed matrix S*T. If full compression
%  is used, it is possible that certain flux modes cannot be decompressed;
%  in this case the corresponding colums in fm_bin have a value of -128.
%  Finally, statistics of the computation are given in the stats table.
%
% [fm_bin,S_unc,id_unc,T,stats] = CALCULATE_FLUX_MODES('.../current_run.mat')
%  In case the automatic subdivsion strategy is used, the intermediate
%  results will be saved after completion of the flux mode computation of
%  each subnetwork ('fm_bin.dat', 'S_unc.dat' and 'current_run.mat' file).
%  If the next iteration gets interrupted, it can be resumed by running
%  CALCULATE_FLUX_MODES with as argument the location of the
%  'current_run.mat' file.
%
%  In order for the out-of-memory detection routine to work optimally, it
%  is strongly advised that virtual memory is turned off by setting the
%  page file size to zero (Windows) or disabling swapping (Linux).
%
%  CALCULATE_FLUX_MODES accepts the following input arguments:
%   'ReactionsUnconstrained' - Indices of the reactions that have to remain
%                              unconstrained. For elementary flux modes no
%                              reactions are unconstrained (default), for
%                              extreme pathways reversible boundary
%                              reactions are unconstrained, for minimal
%                              generating sets all reversible reactions are
%                              unconstrained.
%   'CompoundsUnconstrained' - Indices of the (hub) compounds on which the
%                              steady state constraint will not be imposed.
%   'ReactionOrder' / - Numeric vectors for determining the order in which
%   'CompoundOrder'     reaction/compound constraints are imposed. For
%                       the nullspace approach compound constraints are
%                       resolved first, so CompoundOrder must be a zero
%                       vector and ReactionOrder an all-one vector
%                       (default). For the canonical basis approach
%                       reaction constraints are resolved first, so
%                       ReactionOrder must be all zeros and CompoundOrder
%                       all ones. Custom strategies for determining the
%                       iteration order are possible by giving ReactionOrder
%                       and CompoundOrder any given range of values.
%   'MaxIntermediates' - Maximum number of intermediate flux modes: if
%                        this number is reached the algorithm will stop and
%                        return the intermediate results. Default value is
%                        Inf. In the case an array with two numbers is
%                        given, the algorithm will automatically subdivide
%                        the network for the remaining constraints once the
%                        first threshold MaxIntermediates(1) has been
%                        reached, or if not enough memory is available to
%                        proceed. Subsequently flux modes will be
%                        determined in the subdivided networks and results
%                        will be written to the harddrive (see 'Dir'
%                        option). If the second threshold MaxIntermediates(2)
%                        has been surpassed, the algorithm will stop.
%   'Compression' - Perform flux mode analysis on the compressed
%                   stoichiometric matrix; options are 'full' (default),
%                   'partial' and 'off'. Full compression maximally reduces
%                   the number of reactions, but it may not always be
%                   possible to decompress intermediate flux modes, i.e.
%                   in case the algorithm has stopped because the maximum
%                   number of flux modes has been reached or insufficient
%                   memory was available. In contrast, partial compression
%                   will ensure that also intermediate flux modes can be
%                   decompressed.
%   'StaticOrder' - Strategy that determines the order with which the
%                   constraints are resolved in the double description
%                   method; default 'MaxZero,AbsLexMin'. StaticOrder is
%                   determined before the main iteration. The following
%                   options exist:
%                   'MinZero' - minimum number of zeros
%                   'MaxZero' - maximum number of zeros
%                   'MinNeg' - minimum number of negative coefficients
%                   'MaxNeg' - maximum number of negative coefficients
%                   'MinPos' - minimum number of positive coefficients
%                   'MaxPos' - maximum number of positive coefficients
%                   'MinComb' - minimum number of combinations between
%                               positive and negative coefficients
%                   'MaxComb' - maximum number of combinations between
%                               positive and negative coefficients
%                   'MinMinPosNeg' - minimum number of positive or
%                                    negative coefficients, depending on
%                                    which is the smallest number
%                   'LexMin' - lexicographic order
%                   'LexMax' - inverted lexicographic order
%                   'AbsLexMin' - lexicographic order of the absolute
%                                 coefficient values
%                   'AbsLexMax' - inverted lexicographic order of the
%                                 absolute coefficient values
%                   'Random' - random order
%                   '' - no static order; see 'DynamicOrder'
%                   NB A combination of ordering strategies can be used
%                   by separating each strategy by a comma.
%   'DynamicOrder' - Strategy that determines the order with which the
%                    constraints are resolved in the double description
%                    method; default 'MinComb'. 'DynamicOrder' is
%                    determined during the main iteration of the double
%                    description method and is conditional on 'StaticOrder'.
%                    In case only dynamic ordering has to be used,
%                    'StaticOrder' has to be set to ''. For 'DynamicOrder'
%                    the same options exist as for 'StaticOrder', except
%                    'LexMin', 'LexMax', 'AbsLexMin', 'AbsLexMax' and
%                    'Random'.
%   'Display' - Level of verbosity; options are 'iter' (default) for
%               showing statistics at each iteration, 'on' for showing
%               only general messages, or else 'off' for silent mode.
%   'Tolerance' - Relative threshold under which floating point
%                 coefficients will be rounded off to zero or to the
%                 closest fraction; default value is 1e-10.
%   'Dir' - Directory on the harddrive to store intermediate solutions when
%           using the automatic splitting option or when the uncompressed
%           flux modes do not fit in the memory. In the case multiple media
%           are available with different access speeds, one should pass a
%           cell with two path names of which the first refers to a
%           directory on a fast medium (typically the internal harddrive)
%           and the second to a directory on a slower medium, but with
%           larger storage (e.g. external harddrive).
%   'MinAvailableMemory' - Minimal amount (MBs) of memory that has to be
%                          available for basic operations; default is 1024
%                          MB. Increase this value if out-of-memory
%                          detection does not work and an out-of-memory
%                          error is given.
%   'MaxAvailableMemory' - Maximal total amount (MBs) of memory that is
%                          available for performing flux mode calculation.
%                          Value can be set to Inf if no upper limit of
%                          memory use is required.
%   'MaxThreads' - Maximum number of threads that is used in the adjacency
%                  test. If MaxThreads is larger than the number of logical
%                  processor cores the latter number will be taken instead.
%   'LPsolver' - Name of the linear programming solver that is used for
%                network compression; default 'none'. See the LP_SOLVE
%                function for the options.
%   'SSE4' - Option to indicate whether the current CPU supports the SSE4
%            instruction set; default 'auto', else 'off'. If supported, the
%            pruning function will calculate the number of nonzero elements
%            using the internal popcnt instruction, which is approximately
%            1.5 times faster than lookup tables.
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

tic

if nargin>=2

    param=function_arguments(struct('ReactionsUnconstrained',[],'CompoundsUnconstrained',[],'Display','iter',...
                                    'Compression','full','LPsolver','none','SSE4','auto',...
                                    'StaticOrder','MaxZero,AbsLexMin','DynamicOrder','MinComb',...
                                    'MaxIntermediates',Inf,'Tolerance',1e-10,'ReactionOrder',1,'CompoundOrder',0,...
                                    'Dir','','Condition','on','MinAvailableMemory',1024,'MaxAvailableMemory',Inf,...
                                    'MaxThreads',Inf),varargin);

    if strcmp(lower(param.LPsolver),'none')
        param.LPsolver='';
    end;

    param.StaticOrder=split_str(lower(param.StaticOrder),',');
    for i=1:length(param.StaticOrder)
        if sum(strcmp({'minminposneg','minzero','maxzero','minneg','maxneg','minpos','maxpos','mincomb','maxcomb','lexmin','lexmax','abslexmin','abslexmax','random'},param.StaticOrder{i}))==0
            error(['Unknown option for ''StaticOrder'': ',param.StaticOrder{i}])
        end;
    end;
    param.DynamicOrder=split_str(lower(param.DynamicOrder),',');
    for i=1:length(param.DynamicOrder)
        if sum(strcmp({'minminposneg','minzero','maxzero','minneg','maxneg','minpos','maxpos','mincomb','maxcomb'},param.DynamicOrder{i}))==0
            error(['Unknown option for ''DynamicOrder'': ',param.DynamicOrder{i}])
        end;
    end;

    if sum(strcmp(lower(param.SSE4),{'auto','off'}))>0
        param.SSE4=strcmp(lower(param.SSE4),'auto');
    else
        error(['Unknown option for ''SSE4'': ',param.SSE4,' . . .'])
    end;

    param.MaxThreads=uint16(min(max(round(param.MaxThreads),1),65535));

    param.CompoundCoefficient=ones(size(S,1),1);
    param.ReactionCoefficient=ones(1,size(S,2));



    S_in=S;
    clear S

    if size(S_in,2)~=length(c)
        error('Length of c must be equal to column size of S . . .')
    end;

    if isempty(S_in)
        fm_bin=zeros(size(S_in,2),0,'int8');
        S_unc=zeros(0,0);
        id_unc=zeros(2,0);
        T=eye(size(S_in,2));
        stats=struct('max_coefficient',[],'time',[],'number_FMs',[],'number_constraints',[],'number_candidate_FMs',zeros(3,0),'id',zeros(2,0),'relative_coefficient_distribution',zeros(18,0),'warnings',{cell(0,1)});
        return
    end;

    if 1==1
        % try to make columns in the stoichiometric matrix all integer
        [x,l]=round_fraction(S_in,1e-12);
        S_in(:,l<1e6)=round(S_in(:,l<1e6).*l(ones(size(S_in,1),1),l<1e6));
        param.ReactionCoefficient(:,l<1e6)=param.ReactionCoefficient(:,l<1e6).*l(:,l<1e6);
        x=find(sum(S_in~=round(S_in),1)>0);
        [a,l]=round_fraction(S_in(:,x),1e-6);
        x=x(l<1e6);
        l=l(:,l<1e6);
        for i=1:length(x)
            a=abs(round(S_in(:,x(i))*l(i))/l(i)-S_in(:,x(i)))./abs(S_in(:,x(i)));
            a(isnan(a))=0;
            if max(a)<param.Tolerance
                S_in(:,x(i))=round(S_in(:,x(i))*l(i));
                param.ReactionCoefficient(:,x(i))=param.ReactionCoefficient(:,x(i)).*l(i);
            end;
        end;
    end;

    if strcmp(param.Condition,'on')
        if sum(sum(S_in~=round(S_in)))>0
            % condition stoichiometric matrix
            [S_in,s_left,s_right]=condition_matrix(S_in,'Round2power','on');
            param.CompoundCoefficient=param.CompoundCoefficient.*s_left(:);
            param.ReactionCoefficient=param.ReactionCoefficient.*(s_right(:)');
            clear s_left s_right
        end;
    elseif ~strcmp(param.Condition,'off')
        error(['Unknown option for ''Condition'': ',param.Condition,' . . .'])
    end;

    s=abs(full(S_in(S_in~=0)));
    s=max(max(s))/min(min(s));
    if (s>1e6)&(sum(sum(S_in~=round(S_in)))>0)
        disp(['Warning: large difference between largest and smallest non-zero element in S (|max|/|min| = ',num2str(s,'%0.5g'),')'])
    end;

    if nargin<2
        c=zeros(1,size(S_in,2));
    end;
    if islogical(param.ReactionsUnconstrained)
        param.ReactionsUnconstrained=find(param.ReactionsUnconstrained);
    end;
    if islogical(param.CompoundsUnconstrained)
        param.CompoundsUnconstrained=find(param.CompoundsUnconstrained);
    end;

    param.CompoundOrder=param.CompoundOrder(:);
    param.ReactionOrder=param.ReactionOrder(:);
    if length(param.CompoundOrder)==1
        param.CompoundOrder=ones(size(S_in,1),1)*param.CompoundOrder;
    end;
    if length(param.ReactionOrder)==1
        param.ReactionOrder=ones(size(S_in,2),1)*param.ReactionOrder;
    end;

    if (min([min(param.CompoundsUnconstrained),1])<1)|(max([max(param.CompoundsUnconstrained),1])>size(S_in,1))
        error(['CompoundsUnconstrained indices out of range [1..',num2str(size(S_in,1)),']'])
    end;
    if (min([min(param.ReactionsUnconstrained),1])<1)|(max([max(param.ReactionsUnconstrained),1])>size(S_in,2))
        error(['ReactionsUnconstrained indices out of range [1..',num2str(size(S_in,2)),']'])
    end;
    if sum(param.CompoundsUnconstrained~=round(param.CompoundsUnconstrained))>0
        error(['CompoundsUnconstrained must be integers or logicals'])
    end;
    if sum(param.ReactionsUnconstrained~=round(param.ReactionsUnconstrained))>0
        error(['ReactionsUnconstrained must be integers or logicals'])
    end;
    if size(S_in,2)~=length(c)
        error('Length of c must be equal to the number of columns in S')
    end;
    if size(S_in,1)~=length(param.CompoundOrder)
        error('Length of CompoundOrder must be equal to the number of rows in S')
    end;
    if size(S_in,2)~=length(param.ReactionOrder)
        error('Length of ReactionOrder must be equal to the number of columns in S')
    end;

    if isempty(param.Dir)
        param.Dir=pwd;
    end;
    if ischar(param.Dir)
        param.Dir={param.Dir};
    end;
    if iscell(param.Dir)&(length(param.Dir)==1)
        param.Dir=[param.Dir,param.Dir];
    end;
    if exist(param.Dir{1},'dir')==0
        error(['Invalid directory name: "',param.Dir{1},'" . . .'])
    end;
    if (exist(param.Dir{2},'dir')==0)&(strcmp(param.Dir{1},param.Dir{2})==0)
        error(['Invalid directory name: "',param.Dir{2},'" . . .'])
    end;


    param.Display=lower(param.Display);
    if sum(strcmp(param.Display,{'off','on','iter'}))==0
        error(['Unknown display option: "',param.Display,'" . . .'])
    end;
    param.Compression=lower(param.Compression);
    if sum(strcmp(param.Compression,{'off','partial','full'}))==0
        error(['Unknown compression option: "',param.Compression,'" . . .'])
    end;



    % discretize reaction and compound order
    x=[param.CompoundOrder;param.ReactionOrder];
    x=[x,[1:length(x)]'];
    x=sortrows(x,1);
    x(:,1)=cumsum([1;diff(x(:,1))>1e-10])-1;
    x=sortrows(x,2);
    param.CompoundOrder=x(1:length(param.CompoundOrder),1);
    param.ReactionOrder=x([1:length(param.ReactionOrder)]+length(param.CompoundOrder),1);



    param.Constraints=double(c(:)');
    clear c x





    S=S_in;
    % Calculate reduced row echelon form of S in three phases:
    %  1. rref with hub reactions as pivots; eliminate dependent hub reactions
    %  2. rref with reactions that can be eliminated using compression
    %  3. full rref of S for which CompoundOrder==0

    ReactionOrder=param.ReactionOrder;
    ReactionOrder(param.ReactionsUnconstrained)=Inf;
    CompoundOrder=param.CompoundOrder(:);
    CompoundOrder(param.CompoundsUnconstrained)=Inf;

    % define compound order i_c and reaction order i_r
    [a,i_c]=sort(CompoundOrder);
    i_c=i_c(:);
    S=S(i_c,:);
    i_r=[1:size(S,2)]';


    % phase 1
    j=isinf(ReactionOrder(i_r));
    i_r=[flipud(i_r(j,:));i_r(~j,:)];
    n=sum(j);
    i_CompoundOrder=min(CompoundOrder);
    p=zeros(1,size(S,2))==1;
    while (n>sum(p))&(~isempty(i_CompoundOrder))&(~isinf(i_CompoundOrder))
        if sum(sum(S~=round(S)))==0
            [S,p,i]=rref_fast(S,'ColumnOrder',i_r(1:n),'ExcludeRows',find(CompoundOrder(i_c)>i_CompoundOrder),'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance,'Integer','on');
        else
            [S,p,i]=rref_fast(S,'ColumnOrder',i_r(1:n),'ExcludeRows',find(CompoundOrder(i_c)>i_CompoundOrder),'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance);
        end;
        i_c=i_c(i,:);
        CompoundOrder(i_c(1:sum(p)))=0;
        i_r=i_r([find(p(:,i_r)==1),find(p(:,i_r)==0)]);
        i_CompoundOrder=min(CompoundOrder(CompoundOrder>i_CompoundOrder));
    end;

    % update compound and reaction order and unconstrained indices
    param.ReactionOrder(p)=Inf;
    i=find(~isinf(param.ReactionOrder(param.ReactionsUnconstrained)));
    if (~isempty(i))&(~strcmp(param.Display,'off'))
        disp(['The following reactions could not be left unconstrained: ',num2str(reshape(param.ReactionsUnconstrained(i),1,length(i)))])
    end;
    param.ReactionsUnconstrained(i)=[];
    param.CompoundOrder=CompoundOrder;
    clear CompoundOrder ReactionOrder

    pivot_hub=p;




    % phase 2: compress stoichiometric matrix
    param.T=sparse(1:size(S,2),1:size(S,2),param.ReactionCoefficient);
    if ~strcmp(param.Compression,'off')
        i_excl=ones(1,size(S,2))==0;

        % do not compress zero reactions
        i_excl(full(sparse(1,sort([find(param.Constraints(:)==-1)]),1,1,size(S,2)))~=0)=1==1;

        % do not compress reactions involving more than 15 compounds
        i_excl(sum(S~=0,1)>15)=1==1;

        % do not compress reactions with stoichiometric coefficients > 50
        i_excl(max(abs(S),[],1)>50)=1==1;

        % do not compress reactions with non-integer coefficients
        %i_excl(sum(S~=round(S),1)==0)=[];

        [a,i]=sortrows([sum(S~=round(S),1)'>0,sum(S~=0,1)']);
        i(i_excl(i))=[];

        if sum(sum(S~=round(S)))==0
            [S,p,i]=rref_fast(S,'ColumnOrder',i,'ExcludeRows',find(isinf(param.CompoundOrder(i_c))),'PivotSelection',[param.Compression,' compression'],'NonNegativityConstraints',param.Constraints~=1,'Pivots',pivot_hub,'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance,'Integer','on');
        else
            [S,p,i]=rref_fast(S,'ColumnOrder',i,'ExcludeRows',find(isinf(param.CompoundOrder(i_c))),'PivotSelection',[param.Compression,' compression'],'NonNegativityConstraints',param.Constraints~=1,'Pivots',pivot_hub,'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance);
        end;
        i_c=i_c(i,:);
        pivot_compr=p;
        pivot_compr(pivot_hub==1)=0;

        if ~isempty(param.LPsolver)
            if strcmp(param.Compression,'full')
                [C_blocked,rev_change,C,p]=LP_compression(S((sum(pivot_compr+pivot_hub)+1):end,:),param.Constraints,i_excl,param.Tolerance,param.LPsolver);
                C=[C_blocked;C];
                p(sum(C_blocked,1)>0)=1==1;
            else
                [C,rev_change]=LP_compression(S((sum(pivot_compr+pivot_hub)+1):end,:),param.Constraints,i_excl,param.Tolerance,param.LPsolver);
                p=sum(C,1)>0;
            end;

            if sum(p)>0
                S=[S(1:sum(pivot_compr+pivot_hub),:);C;S((sum(pivot_compr+pivot_hub)+1):end,:)];
                i_c=[i_c(1:sum(pivot_compr+pivot_hub),:);size(S_in,1)+[1:size(C,1)]';i_c((sum(pivot_compr+pivot_hub)+1):end,:)];
                S_in=[S_in;C];
                param.T(:,rev_change==-1)=-param.T(:,rev_change==-1);
                S(:,rev_change==-1)=-S(:,rev_change==-1);
                S_in(:,rev_change==-1)=-S_in(:,rev_change==-1);
                param.Constraints(rev_change~=0)=0;

                param.CompoundCoefficient=[param.CompoundCoefficient;ones(size(C,1),1)];
                param.CompoundOrder=[param.CompoundOrder;zeros(size(C,1),1)];
                if sum(sum(S~=round(S)))==0
                    [S,p,i]=rref_fast(S,'ColumnOrder',[find(pivot_hub),find(pivot_compr),find(p)],'ExcludeRows',size(C,1)+sum(pivot_compr+pivot_hub)+[1:(size(S,1)-size(C,1)-sum(pivot_compr+pivot_hub))],'Pivots',pivot_compr|pivot_hub,'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance,'Integer','on');
                else
                    [S,p,i]=rref_fast(S,'ColumnOrder',[find(pivot_hub),find(pivot_compr),find(p)],'ExcludeRows',size(C,1)+sum(pivot_compr+pivot_hub)+[1:(size(S,1)-size(C,1)-sum(pivot_compr+pivot_hub))],'Pivots',pivot_compr|pivot_hub,'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance);
                end;
                [a,i]=sort(i);
                S=S(i,:);
                pivot_compr(p)=1;
                pivot_compr(pivot_hub==1)=0;
            elseif sum(abs(rev_change))>0
                param.T(:,rev_change==-1)=-param.T(:,rev_change==-1);
                S(:,rev_change==-1)=-S(:,rev_change==-1);
                S_in(:,rev_change==-1)=-S_in(:,rev_change==-1);
                param.Constraints(rev_change~=0)=0;
            end;
        end;

        T=sparse(null_rref(S,pivot_compr));

        param.T=param.T*T;
        param.CompoundOrder(i_c(1:sum(p)))=0;

        if ~strcmp(param.Display,'off')
            if size(T,1)>size(T,2)
                disp(['Compression on: reduced the number of reactions from ',num2str(size(T,1)),' to ',num2str(size(T,2))])
            else
                disp('Compression on: the number of reactions could not be reduced  ')
            end;
        end;
    else
        T=sparse(1:size(S,2),1:size(S,2),1);
        pivot_compr=zeros(1,size(S,2));
        if ~strcmp(param.Display,'off')
            disp('Compression off')
        end;
    end;



    % phase 3
    % calculate rref for all rows in S with order index 0 and save result for determining the initial R


    [a,i]=sortrows([[1:size(S,1)]'>sum(pivot_hub+pivot_compr),param.CompoundOrder(i_c,:)]);
    S=S(i,:);
    i_c=i_c(i,:);
    S1=S((sum(pivot_hub+pivot_compr)+1):end,:);
    a=max(abs(S1'),[],2)./min(abs(S1').*(S1'~=0)./(S1'~=0),[],2);
    a(isnan(a))=0;
    [x,i_r]=sortrows([-param.ReactionOrder(:),a>=4,full(sum(S1~=0,1)')>=10,full(sum(S1~=0,1)')>=5,param.Constraints(:)~=1,full(sum(S1~=0,1)')]);
    clear S1

    if sum(sum(S~=round(S)))==0
        [S,p,i,j]=rref_fast(S,'ColumnOrder',i_r,'ExcludeRows',find((param.CompoundOrder(i_c,:)~=0)&(sum(S(:,pivot_hub|pivot_compr)~=0,2)==0)),'Pivots',pivot_hub+pivot_compr,'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance,'Integer','on');
    else
        [S,p,i,j]=rref_fast(S,'ColumnOrder',i_r,'ExcludeRows',find((param.CompoundOrder(i_c,:)~=0)&(sum(S(:,pivot_hub|pivot_compr)~=0,2)==0)),'Pivots',pivot_hub+pivot_compr,'ReferenceMatrix',S_in(i_c,:),'Tolerance',param.Tolerance);
    end;
    i_c=i_c(i,:);
    i_r=i_r([find(p(:,i_r)==1),find(p(:,i_r)==0)]);
    pivot_DD=p;
    pivot_DD((pivot_hub+pivot_compr)>0)=0;
    if ~isempty(j)
        S((j(:)==0)&(~isinf(param.CompoundOrder(i_c,:))),:)=0;
    end;

    clear S_in

    % eliminate compressed part of S
    i=sparse(T~=0);
    i(pivot_compr~=0,:)=0;
    S1=S*i;
    S1(sum(S(:,pivot_compr~=0)~=0,2)==1,:)=0;
    S(sum(S(:,pivot_compr~=0)~=0,2)==1,:)=0;
    S(:,pivot_compr~=0)=[];
    if max(max(abs(S-S1)))>param.Tolerance
        error('Invalid compression matrix T . . .')
    end;

    % Edit #1
    disp('Compressed S:');
    disp(S);
    clear S1

    i=zeros(sum(sum(T~=0)),3);
    [i(:,1),i(:,2),i(:,3)]=find(T);

    % reassign reaction specific parameters
    a=sortrows([i(:,2),param.Constraints(:,i(:,1))'],[1,(1:size(param.Constraints,1)+1)]);
    a(find(diff(a(:,1))==0)+1,:)=[];
    param.Constraints=a(:,2:end)';
    i=sparse(T~=0);
    i(pivot_compr~=0,:)=0;

    pivot=full((pivot_hub+pivot_DD)*i)==1;
    clear pivot_hub pivot_DD pivot_compr

    param.ReactionOrder=i'*param.ReactionOrder;
    n=length(param.ReactionsUnconstrained);
    [i1,i2,i3]=find(i'*sparse(param.ReactionsUnconstrained,1:n,1,size(T,1),n));
    i=sortrows([i2(:),i1(:)]);
    param.ReactionsUnconstrained=i(:,2);
    clear i1 i2 i3 n i j
    a=isinf(param.ReactionOrder);
    if sum(~a(param.ReactionsUnconstrained))>0
        error('Invalid compression matrix T . . .')
    end;
    a(param.ReactionsUnconstrained)=[];
    if sum(a)>0
        error('Invalid compression matrix T . . .')
    end;



    % expand reversible reactions that are not indicated as hub reactions
    irrev2rev_mapping=(param.Constraints(:)'==1)&(~isinf(param.ReactionOrder(:)'));
    S=[S,-S(:,irrev2rev_mapping)];
    pivot=[pivot,zeros(1,sum(irrev2rev_mapping))==1];
    param.Constraints=[param.Constraints,zeros(1,sum(irrev2rev_mapping))+NaN];
    param.ReactionOrder=param.ReactionOrder(:)';
    param.ReactionOrder=[param.ReactionOrder,param.ReactionOrder(:,irrev2rev_mapping)];

    param.S_red=S;
    param.p_red=pivot;
    param.i_r=i_r;
    param.i_c=i_c;
    clear S pivot i_r i_c



elseif (nargin==1)&(ischar(S))
    % resume network splitting session

    if exist(S,'file')
        load(S)
        param.Dir=fileparts(S);
        if isempty(param.Dir)
            param.Dir=pwd;
        end;
        param.Dir={param.Dir,param.Dir};
        fm_bin=[param.Dir{2},'/fm_bin.dat'];
        S_unc=[param.Dir{2},'/S_unc.dat'];
    else
        error(['File does not exist: ',S,' . . .'])
    end;

else
    error('Expected MAT file to resume network splitting session . . . ')
end;

% disp('****Now Here!<2>****')
% pause

if length(param.MaxIntermediates)==1
    % disp('****Now Here!<3>****')
    % pause

    [param.S_red,param.p_red,param.i_c]=update_S(param.S_red,param.p_red,irrev2rev_mapping,[],[],param.i_c,param.Tolerance,'',(sum(param.S_red(:,isinf(param.ReactionOrder))~=0,2)>0)|isinf(param.CompoundCoefficient));
    S_th=param.S_red*0;
    S_th(param.i_c,:)=param.S_red;

    [A,i_A,i_A_resolved,order_A,N,param]=determine_initial_solution(param.S_red,param.p_red,param.i_c,param.i_r,param);


    exit_str='';
    initialise_DD(A,i_A,order_A,N,S_th,param,i_A_resolved);
    exit_str=main_DD;
    if ~strcmp(param.Display,'off')
        if ~isempty(exit_str)
            disp(exit_str)
        else
            disp('Double description method finished successfully')
        end;
    end;

    % disp('****Now Here!<4>****')
    % pause

    [fm_bin,S_unc,id_unc,stats]=get_DD_results(irrev2rev_mapping==1);
    i=find(irrev2rev_mapping);
    j=(stats.id(2,:)>length(irrev2rev_mapping))&(stats.id(1,:)==2);
    stats.id(2,j)=-i(stats.id(2,j)-length(irrev2rev_mapping));
    j=(id_unc(2,:)>length(irrev2rev_mapping))&(id_unc(1,:)==2);
    id_unc(2,j)=-i(id_unc(2,j)-length(irrev2rev_mapping));

    T=sparse(param.T);
    i=find((id_unc(1,:)==2)&(id_unc(2,:)<0));
    if ~isempty(i)
        id_unc(2,i)=-id_unc(2,i);
        T(:,i)=-T(:,i);
    end;

    disp(['Number of flux modes: ',num2str(stats.number_FMs(end))])
    disp(['Total computation time: ',num2str(stats.time(end)),' s'])

    if ischar(fm_bin)&ischar(S_unc)
        if exist(fm_bin,'file')&exist(S_unc,'file')
            % if flux modes are written to harddrive, write all results to harddrive
            if exist([param.Dir{2},'/fm.mat'],'file')
                delete([param.Dir{2},'/fm.mat']);
            end;
            save([param.Dir{2},'/fm.mat'],'id_unc','T','stats','-V6');
            id_unc=[param.Dir{2},'/fm.mat'];
            T=[];
            stats=[];
        else
            error('???')
        end;
    end;

elseif length(param.MaxIntermediates)==2

    if nargin>=2
        i_iter=0;
        is_zero=zeros(size(param.S_red,2),1)==1;
        is_nonzero=zeros(size(param.S_red,2),1)==1;
        n=0;
        n1=0;
        exit_str='';
    end;

    if max(param.CompoundOrder)>min(param.ReactionOrder)
        error('Nullspace approach must be used for splitting . . .')
    end;

    while (~isempty(is_zero))&(isempty(exit_str))
        i_iter=i_iter+1;
        i_stack=size(is_zero,2);

        if strcmp(param.Display,'iter')
            disp(['Splitting strategy iteration ',num2str(i_iter),';  stack size ',num2str(size(is_zero,2)),';  # intermediate flux modes: ',num2str(n)])
        end;
        param_local=param;
        param_local.MaxIntermediates=param_local.MaxIntermediates(1);

        i=update_S(param_local.S_red,param_local.p_red,irrev2rev_mapping,is_zero(:,i_stack)==1,is_nonzero(:,i_stack)==1,param_local.i_c,param_local.Tolerance,param_local.LPsolver,(sum(param.S_red(:,isinf(param.ReactionOrder))~=0,2)>0)|isinf(param.CompoundCoefficient));
        is_zero(i,i_stack)=1;

        param_local.Constraints(is_zero(:,i_stack)==1)=-1;
        S_th=param_local.S_red*0;
        S_th(param_local.i_c,:)=param_local.S_red;

        [A,i_A,i_A_resolved,order_A,N,param_local]=determine_initial_solution(param_local.S_red,param_local.p_red,param_local.i_c,param_local.i_r,param_local);


        i=find(i_A(1,:)==1);
        order_A(1,i(is_zero(i_A(2,i),i_stack)==1))=min(order_A(1,:))-1;
        order_A(1,i(is_nonzero(i_A(2,i),i_stack)==1))=Inf;
        i=find(i_A(1,:)==2);
        order_A(1,i(is_zero(i_A(2,i),i_stack)==1))=min(order_A(1,:))-1;
        order_A(1,i(is_nonzero(i_A(2,i),i_stack)==1))=Inf;

        exit_str1='';
        initialise_DD(A,i_A,order_A,N,S_th,param_local,i_A_resolved);
        exit_str1=main_DD;
        if ~strcmp(param.Display,'off')
            if ~isempty(exit_str1)
                disp(exit_str1)
            else
                disp('Double description method finished successfully')
            end;
        end;
        [x,id_unc_local]=filter_splitting_results(S_th);

        j=zeros(1,size(id_unc_local,2))==1;
        i=full(sparse(param.CompoundsUnconstrained,1,1,size(S_th,1),1))~=0;
        j(id_unc_local(1,:)==1)=i(id_unc_local(2,id_unc_local(1,:)==1));
        i=full(sparse(param.ReactionsUnconstrained,1,1,size(S_th,2),1))~=0;
        j(id_unc_local(1,:)==2)=i(id_unc_local(2,id_unc_local(1,:)==2));


        if i_iter>1
            % check whether id_S_unc has not changed
            id_unc1=id_unc_local(:,j);
            if size(id_unc1,2)~=size(id_unc,2)
                error('????')
            elseif sum(sum(id_unc~=id_unc1))>0
                error('????')
            end;
        else
            id_unc=id_unc_local(:,j);

            fm_bin=[param.Dir{2},'/fm_bin.dat'];
            S_unc=[param.Dir{2},'/S_unc.dat'];
            if exist(fm_bin,'file')
                delete(fm_bin);
            end;
            if exist(S_unc,'file')
                delete(S_unc);
            end;
            if exist([param.Dir{2},'/fm.mat'],'file')
                delete([param.Dir{2},'/fm.mat']);
            end;
        end;
        [fm_bin_local,S_unc_local,id_unc_local,stats_local]=get_DD_results(irrev2rev_mapping==1,1==1,find(~j));

        i=find(irrev2rev_mapping);
        ii=(stats_local.id(2,:)>length(irrev2rev_mapping))&(stats_local.id(1,:)==2);
        stats_local.id(2,ii)=-i(stats_local.id(2,ii)-length(irrev2rev_mapping));

        n1=n1+x;
        n=n+stats_local.number_FMs(end);
        clear fm_bin_local
        clear S_unc_local
        if (n>param.MaxIntermediates(2))
            exit_str='Maximum number of intermediate flux modes in total network has been reached ';
        end;

        stats_local.stack_size=size(is_zero,2);
        stats_local.number_FMs_total=[n;n1];
        if i_iter>1
            stats(i_iter)=stats_local;
        else
            stats=stats_local;
        end;

        if sum(id_unc_local(1,~j)==1)>0
            error('???')
        end;

        i=id_unc_local(2,~j);
        i(is_nonzero(i,i_stack))=[];
        if isempty(N)
            i=[];
        end;
        if (sum(is_zero(i,i_stack))>0)
            exit_str='Could not resolve all zero constraints ';
        end;
        is_zero_new=is_zero(:,ones(1,length(i))*i_stack);
        is_zero_new(i,:)=eye(length(i))==1;
        is_nonzero_new=is_nonzero(:,ones(1,length(i))*i_stack);
        is_nonzero_new(i,:)=triu(ones(length(i),length(i)),1)==1;



        is_zero(:,i_stack)=[];
        is_nonzero(:,i_stack)=[];
        %is_zero=[is_zero,is_zero_new];
        %is_nonzero=[is_nonzero,is_nonzero_new];
        is_zero=[is_zero_new,is_zero];
        is_nonzero=[is_nonzero_new,is_nonzero];
        disp(' ')

        save([param.Dir{2},'/current_run.mat'],'i_iter','is_zero','is_nonzero','n','n1','exit_str','param','irrev2rev_mapping','stats','id_unc','-V6');
    end;

    if ~strcmp(param.Display,'off')
        if ~isempty(exit_str)
            disp(exit_str)
        else
            disp('Automatic splitting method finished successfully')
        end;
    end;

    T=sparse(param.T);
    i=find((id_unc(1,:)==2)&(id_unc(2,:)<0));
    if ~isempty(i)
        id_unc(2,i)=-id_unc(2,i);
        T(:,i)=-T(:,i);
    end;

    save([param.Dir{2},'/fm.mat'],'id_unc','T','stats','-V6');
    disp(['Number of flux modes: ',num2str(stats(end).number_FMs_total(1,end))])
    disp(['Total computation time: ',num2str(stats(end).time(end)),' s'])
    id_unc=[param.Dir{2},'/fm.mat'];
    T=[];
    stats=[];
end;







function [S,p,i_c,i_zero_out]=update_S(S_in,p_in,irrev2rev_mapping,i_zero,i_nonzero,i_c_in,tol,LPsolver,row_excl)

S=S_in;
i_c=i_c_in;
p=p_in~=0;

if islogical(i_zero)
    i_zero=find(i_zero);
end;
if islogical(i_nonzero)
    i_nonzero=find(i_nonzero);
end;

% determine initial C
j=full([sparse(1,find(irrev2rev_mapping),length(irrev2rev_mapping)+[1:sum(irrev2rev_mapping)],1,length(irrev2rev_mapping)),find(irrev2rev_mapping(:)')]);
i=j(i_nonzero);
i(i==0)=[];
i=sort([i(:);i_zero(:)]);
i(diff(i)==0)=[];
C=full(sparse(1:length(i),i,1,length(i),size(S,2)));

% extend C by global constraints
if ~isempty(LPsolver)
    % only use LP compression if a stable solver is used, in order to prevent cumulative memory leakage

    i=find(sum(C~=0,1)>0);

    S1=S(:,1:length(irrev2rev_mapping));
    rev1=double(irrev2rev_mapping(:)');
    k=[1:size(S1,2)];

    S1(:,i(j(i)>length(irrev2rev_mapping)))=-S1(:,i(j(i)>length(irrev2rev_mapping)));
    rev1(i(j(i)==0))=-1;
    rev1(j(i(j(i)~=0)))=0;
    rev1(i(j(i)~=0))=0;
    rev1=rev1(:,1:length(irrev2rev_mapping));

    [C1,rev_change]=LP_compression(S1(~row_excl,:),rev1,[],tol,LPsolver);


    C1=find(sum(C1,1)>0);
    C1=[C1(:)',j(:,C1(j(C1)~=0))];
    C1=sort([C1,j(:,rev_change==1),find(rev_change==-1)]);
    C1(C1==0)=[];
    C1(diff(C1)==0)=[];
    C1(sum(C(:,C1)~=0,1)>0)=[];
    C=[C;full(sparse(1:length(C1),C1,1,length(C1),size(C,2)))];
end;



if nargout==1
    if sum(sum(C(:,i_nonzero)~=0))==0
        S=sum(C~=0,1)>0;
    else
        S=ones(1,size(S,2))==1;
    end;
    return
end;

% update S i_c
if ~isempty(C)
    if sum(sum(C(:,i_nonzero)~=0))>0
        S=eye(size(S,2));
        p=ones(1,size(S,2))==1;
        i_c=[1:size(S,2)]';
        i_zero_out=ones(1,size(S,2))==1;
    else
        i_c=[max(i_c)+[1:size(C,1)]';i_c];
        i_excl=[zeros(size(C,1),1)==1;sum(S(:,p)~=0,2)==0];
        S(:,sum(C~=0,1)>0)=0;
        S=[C;S];
        p(sum(C~=0,1)>0)=1==1;

        s=ones(size(S,2),1);
        s(i_nonzero)=0;
        [a,i]=sortrows([s,sum(S~=0,1)',sum(abs(S),1)']);
        if sum(sum(S~=round(S)))==0
            [S,p1,i]=rref_fast(S,'ColumnOrder',i,'ExcludeRows',find(i_excl|(sum(S(:,p)~=0,2)>0)),'Tolerance',tol,'Integer','on');
        else
            [S,p1,i]=rref_fast(S,'ColumnOrder',i,'ExcludeRows',find(i_excl|(sum(S(:,p)~=0,2)>0)),'Tolerance',tol);
        end;
        [a,i]=sort(i);
        S=S(i,:);
        p(p1)=1==1;
    end;
end;
i_zero_out=sum(C~=0,1)>0;




function [A,i_A,i_A_resolved,order_A,N,param]=determine_initial_solution(S,pivot,i_c,i_r,param_in)

param=param_in;


if sum(pivot)>0
    % calculate null space matrix
    if 1==1
        % do not enforce lcm in case of >1 pivot coefficients
        i=sum(S(:,pivot==1)~=0,2)>0;
        N=zeros(size(S,2),sum(pivot==0));
        N(pivot==0,:)=eye(sum(pivot==0));
        N(pivot==1,:)=-(S(i,pivot==1)~=0)'*S(i,pivot==0);
        s=ones(1,size(S,2));
        s(:,pivot==1)=sum(S(i,pivot==1),1);
        if max(max(abs(S(i,:)*diag(1./s)*N)./(abs(S(i,:))*diag(1./s)*abs(N))))>param.Tolerance
            error('Nullspace matrix could not be calculated . . .')
        end;
        param.T=param.T*sparse(1:size(param.T,2),1:size(param.T,2),1./s(:,1:size(param.T,2)));
    else
        N=full(null_rref(S,pivot));
    end;

    N(:,sum(N((param.Constraints(:)==-1)&(~pivot(:)),:)~=0,1)>0)=[]

    i_A_resolved=i_c(sum(S(:,pivot)~=0,2)==1);
    i_A_resolved=[ones(1,length(i_A_resolved));i_A_resolved(:)'];
else
    N=sparse(1:size(S,2),1:size(S,2),1);
    i_A_resolved=zeros(2,0);
end;


i_r=[1:size(S,2)];
i_r(param.ReactionsUnconstrained)=[];
i_A_resolved=[i_A_resolved,[2*ones(1,sum(~pivot(i_r)));i_r(:,~pivot(i_r))]];
i_r(~pivot(i_r))=[];

i=full(sparse(param.CompoundsUnconstrained,1,1,size(S,1),1))~=0;
i=i(i_c,:);
j=full(sum(S~=0,2))~=0;
j(sum(S(:,pivot),2)~=0)=1==0;
A=full([S(j&i,~pivot)*sparse(N(~pivot,:));N(param.ReactionsUnconstrained,:);S(j&(~i),~pivot)*sparse(N(~pivot,:));N(i_r,:)]);


% perform fraction rounding
x=find(sum(A~=round(A),1)>0);
[a,l]=round_fraction(A(:,x),1e-12);
x=x(l<1e6);
l=l(:,l<1e6);
A(:,x)=round(A(:,x).*l(ones(size(A,1),1),:));
A=A';

param.all_integer=sum(sum((A~=round(A))))==0;
if ~strcmp(param.Display,'off')
    if param.all_integer
        disp('Integer mode on')
    else
        disp('Integer mode off')
    end;
end;


i_A=[[ones(1,sum(j&i));i_c(j&i,:)'],[2*ones(1,length(param.ReactionsUnconstrained));param.ReactionsUnconstrained(:)'],...
     [ones(1,sum(j&(~i)));i_c(j&(~i),:)'],[2*ones(1,length(i_r));i_r(:)']];



order_A=zeros(1,size(A,2));
order_A(:,i_A(1,:)==1)=param.CompoundOrder(i_A(2,i_A(1,:)==1));
order_A(:,i_A(1,:)==2)=param.ReactionOrder(i_A(2,i_A(1,:)==2));
x=[];

n_j=[sum(A>0,1);sum(A==0,1);sum(A<0,1)];
[a,j]=sort(sum(A~=0,2));
j(abs(a)==1)=[];
for i=1:length(param.StaticOrder)
    if strcmp(param.StaticOrder{i},'minzero')
        x=n_j(2,:);
    elseif strcmp(param.StaticOrder{i},'maxzero')
        x=-n_j(2,:);
    elseif strcmp(param.StaticOrder{i},'minneg')
        x=n_j(3,:);
    elseif strcmp(param.StaticOrder{i},'maxneg')
        x=-n_j(3,:);
    elseif strcmp(param.StaticOrder{i},'minpos')
        x=n_j(1,:);
    elseif strcmp(param.StaticOrder{i},'maxpos')
        x=-n_j(1,:);
    elseif strcmp(param.StaticOrder{i},'mincomb')
        x=n_j(1,:).*n_j(3,:);
    elseif strcmp(param.StaticOrder{i},'maxcomb')
        x=-n_j(1,:).*n_j(3,:);
    elseif strcmp(param.StaticOrder{i},'minminposneg')
        x=min(n_j([1,3],:),[],1);
    elseif strcmp(param.StaticOrder{i},'lexmin')
        x=A(j,:);
    elseif strcmp(param.StaticOrder{i},'lexmax')
        x=-A(j,:);
    elseif strcmp(param.StaticOrder{i},'abslexmin')
        x=abs(A(j,:));
    elseif strcmp(param.StaticOrder{i},'abslexmax')
        x=-abs(A(j,:));
    elseif strcmp(param.StaticOrder{i},'random')
        x=rand(1,size(A,2));
    else
        error('???')
    end;
    order_A=[order_A;x];
end;

N=full(N)~=0;
N(pivot,:)=0;






function [n,i_A]=filter_splitting_results(S)

global A R i_A param pointer_type

n=size(R,2);

j=zeros(1,size(i_A,2))==1;
i=full(sparse(param.CompoundsUnconstrained,1,1,size(S,1),1))~=0;
j(i_A(1,:)==1)=i(i_A(2,i_A(1,:)==1));
i=full(sparse(param.ReactionsUnconstrained,1,1,size(S,2),1))~=0;
j(i_A(1,:)==2)=i(i_A(2,i_A(1,:)==2));

j=find(~j);


if iscell(A)
    A_size=A{2};
else
    A_size=size(A);
end;
if A_size(1)>0
    buffer_size=get_free_memory;
    if ((buffer_size<A_size(2)+param.MinAvailableMemory*1048576)&(~iscell(R)))
        buffer_size=buffer_size+size(R,1)*size(R,2)*4;
        R=matrix_indexing(R,[param.Dir{1},'/R.dat']);
    end;
    if ((buffer_size<A_size(2)+param.MinAvailableMemory*1048576)&(~iscell(A)))
        buffer_size=buffer_size+size(A,1)*size(A,2)*8;
        A=matrix_indexing(A,[param.Dir{2},'/A.dat']);
    end;
    i=true(1,A_size(2));
    mem_block=round(0.25*buffer_size/8);
    n_block=floor(mem_block/A_size(1));


    if iscell(A)
        fid=fopen(A{1},'r');
    end;
    m=[1-n_block,0];
    while m(2)<A_size(2)
        m=min(m+n_block,A_size(2));
        if iscell(A)
            A_sub=fread(fid,[A_size(1),m(2)-m(1)+1],'double=>double');
        else
            A_sub=A(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));
        end;
        i(feval(pointer_type,m(1)):feval(pointer_type,m(2)))=sum(A_sub<=0,1)==0;
        clear A_sub
    end
    if iscell(A)
        fclose(fid);
        flush_file_cache(A{1});
    end;


    if iscell(R)
        R1=matrix_indexing(R,[param.Dir{1},'/R1.dat'],Inf,i);
        delete(R{1})
        movefile(R1{1},R{1})
        R1{1}=R{1};
        R=R1;
        buffer_size=get_free_memory;
        if buffer_size>2*4*R{2}(1)*R{2}(2)
            a=R{1};
            R=matrix_indexing(R);
            delete(a)
        end;
    else
        buffer_size=get_free_memory;
        if buffer_size>2*4*size(R,1)*sum(i)
            R=matrix_indexing(R,'',Inf,i);
        else
            R1=matrix_indexing(R,[param.Dir{1},'/R.dat'],Inf,i);
            R=R1;
        end;
    end;



    if iscell(A)
        A1=matrix_indexing(A,[param.Dir{2},'/A1.dat'],Inf,i);
        delete(A{1})
        movefile(A1{1},A{1})
        A1{1}=A{1};
        A=A1;
        buffer_size=get_free_memory;
        if buffer_size>2*8*A{2}(1)*A{2}(2)
            a=A{1};
            A=matrix_indexing(A);
            delete(a)
        end;
    else
        buffer_size=get_free_memory;
        if buffer_size>2*8*size(A,1)*sum(i)
            A=matrix_indexing(A,'',Inf,i);
        else
            A1=matrix_indexing(A,[param.Dir{2},'/A.dat'],Inf,i);
            A=A1;
        end;
    end;

end;






function  exit_str=main_DD

global param j pointer_type


% main iteration of the DD method
[n_j,i_current_column,is_zero]=get_next_column;
if ischar(j)
    exit_str=j;
else
    exit_str='';
end;
while (~isempty(i_current_column))&(isempty(exit_str))
    if n_j(1)*n_j(3)>0
        % calculate adjacent j1-j3 combinations
        j_comb=find_adjacent_rays(n_j,n_j(2)+double(is_zero==0)*n_j(1));
    else
        j_comb=zeros(2,0,pointer_type);
    end;

    if ~ischar(j_comb)
        % update DD-pair
        update_DD(n_j,j_comb,i_current_column,is_zero);

        [n_j,i_current_column,is_zero]=get_next_column;
        if ischar(j)
            exit_str=j;
        end;
    else
        exit_str=j_comb;
    end;
end;




function [fm_bin,S_unc,i_S_unc,stats_out]=get_DD_results(irrev2rev_mapping,automatic_splitting,i_A_del);

global A R stats i_A order_A param pointer_type


if nargin<2
    automatic_splitting=1==0;
end;
if nargin<3
    i_A_del=[];
end;
if islogical(i_A_del)
    i_A_del=find(i_A_del);
end;




if iscell(A)
    A_size=A{2};
else
    A_size=size(A);
end;
if iscell(R)
    R_size=R{2};
else
    R_size=size(R);
end;
buffer_size=get_free_memory;
if (buffer_size<(A_size(2)+param.MinAvailableMemory*1048576))&(~iscell(R))
    buffer_size=buffer_size+size(R,1)*size(R,2)*4;
    R=matrix_indexing(R,[param.Dir{1},'/R.dat']);
end;
i_inc=true(1,A_size(2));
buffer_size=buffer_size-A_size(2);

% delete zero cycles in reversible reactions
i_rev=find(irrev2rev_mapping);
i_rev=[i_rev(:),[1:length(i_rev)]'+length(irrev2rev_mapping)];
x=find(i_A(1,:)==2);
x=full(sparse(1,i_A(2,x),x,1,max(max(max([1,1;i_rev])),max([1,i_A(2,x)]))));


n_block=min(floor(0.25*buffer_size/(4*R_size(1)+8*A_size(1))),A_size(2));


if iscell(R)
    fid=fopen(R{1},'r');
end;
if iscell(A)
    fid1=fopen(A{1},'r');
end;
m=[1-n_block,0];
n=0;
while m(2)<R_size(2)
    m=min(m+n_block,R_size(2));

    if iscell(A)
        A_sub=fread(fid1,[A_size(1),m(2)-m(1)+1],'double=>double')';
    else
        A_sub=A(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)))';
    end;
    if iscell(R)
        R_sub=fread(fid,[R_size(1),m(2)-m(1)+1],'uint32=>uint32')';
    else
        R_sub=R(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)))';
    end;

    for i=1:size(i_rev,1)
        if (x(i_rev(i,1))==0)&(x(i_rev(i,2))==0)
            j=floor((i_rev(i,:)-1)/32);
            j=[j+1;i_rev(i,:)-32*j-1];
            j=find((bitand(R_sub(:,j(1,1)),uint32(2^j(2,1)))>0)&(bitand(R_sub(:,j(1,2)),uint32(2^j(2,2)))>0));
        elseif (x(i_rev(i,1))==0)&(x(i_rev(i,2))~=0)
            j=floor((i_rev(i,1)-1)/32);
            j=[j+1;i_rev(i,1)-32*j-1];
            j=find((bitand(R_sub(:,j(1,1)),uint32(2^j(2,1)))>0)&(A_sub(:,x(i_rev(i,2)))~=0));
        elseif (x(i_rev(i,1))~=0)&(x(i_rev(i,2))==0)
            j=floor((i_rev(i,2)-1)/32);
            j=[j+1;i_rev(i,2)-32*j-1];
            j=find((bitand(R_sub(:,j(1,1)),uint32(2^j(2,1)))>0)&(A_sub(:,x(i_rev(i,1)))~=0));
        elseif (x(i_rev(i,1))~=0)&(x(i_rev(i,2))~=0)
            j=find((A_sub(:,x(i_rev(i,1)))~=0)&(A_sub(:,x(i_rev(i,2)))~=0));
        end;
        if length(j)==1
            i_inc(m(1)+j-1)=1==0;
        elseif length(j)>1
            if ~strcmp(param.Display,'off')
                disp('Incorrect EFMs: zero cycle fluxes should only occur once . . .')
            end;
        end;
    end;
    clear A_sub R_sub
end;


if iscell(R)
    fclose(fid);
    flush_file_cache(R{1});
end;
if iscell(A)
    fclose(fid1);
    flush_file_cache(A{1});
end;



if iscell(A)
    A1=A;
    movefile(A1{1},[param.Dir{2},'/A1.dat']);
    A1{1}=[param.Dir{2},'/A1.dat'];
    A=matrix_indexing(A1,[param.Dir{2},'/A.dat'],Inf,i_inc);
    delete(A1{1})
    clear A1
elseif buffer_size<(2*8*A_size(1)*A_size(2)+param.MinAvailableMemory*1048576)
    A_file=matrix_indexing(A,[param.Dir{2},'/A.dat'],Inf,i_inc);
    A=[];
    A=matrix_indexing(A_file);
    delete(A_file{1})
    clear A_file
else
    A=matrix_indexing(A,'',Inf,i_inc);
end;

if iscell(R)
    R1=R;
    movefile(R1{1},[param.Dir{1},'/R1.dat']);
    R1{1}=[param.Dir{1},'/R1.dat'];
    R=matrix_indexing(R1,[param.Dir{1},'/R.dat'],Inf,i_inc);
    delete(R1{1})
    clear R1
elseif buffer_size<(2*4*R_size(1)*R_size(2)+param.MinAvailableMemory*1048576)
    R_file=matrix_indexing(R,[param.Dir{1},'/R.dat'],Inf,i_inc);
    R=[];
    R=matrix_indexing(R_file);
    delete(R_file{1})
    clear R_file
else
    R=matrix_indexing(R,'',Inf,i_inc);
end;
R_size(2)=sum(i_inc);
A_size(2)=sum(i_inc);
clear i_inc


stats.number_FMs=[stats.number_FMs,R_size(2)];
stats.number_constraints=[stats.number_constraints,A_size(1)];

% unpack R
buffer_size=get_free_memory;
n_block=min(floor(0.25*buffer_size/(4*R_size(1) + 8*A_size(1) + 2*size(param.T,1) + 4*2*(strcmp(pointer_type,'double')+1))),A_size(2));

if automatic_splitting|(n_block<R_size(2))
    fm_bin={[param.Dir{2},'/fm_bin.dat'],[size(param.T,1),R_size(2)],'int8'};
    if ~automatic_splitting
        if exist(fm_bin{1},'file')
            delete(fm_bin{1});
        end;
    end;
else
    fm_bin=[];
end;

if iscell(fm_bin)
    fid=fopen(fm_bin{1},'A');
end;
if iscell(R)
    fid1=fopen(R{1},'r');
end;
if iscell(A)
    fid2=fopen(A{1},'r');
end;


if R_size(2)>0
    T=sparse(param.T);
    T=[T,-T(:,irrev2rev_mapping)];
    I=zeros(sum(sum(T~=0)),3);
    [I(:,1),I(:,2),I(:,3)]=find(sign(T));

    m=[1-n_block,0];
    while m(2)<R_size(2)
        m=min(m+n_block,R_size(2));

        fm_bin_sub=ones(m(2)-m(1)+1,size(param.T,1),'int8');
        fm_bin_sub(:)=0;
        if iscell(R)
            R_sub=fread(fid1,[R_size(1),m(2)-m(1)+1],'uint32=>uint32')';
        else
            R_sub=R(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)))';
        end;
        if iscell(A)
            A_sub=fread(fid2,[A_size(1),m(2)-m(1)+1],'double=>double')';
        else
            A_sub=A(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)))';
        end;

        i_inc=zeros(size(param.T,1),1);
        i_inc1=zeros(size(param.T,1),1);
        i_err1=false(m(2)-m(1)+1,1);

        for i=1:size(I,1)
            ii=find((I(i,2)==i_A(2,:))&(i_A(1,:)==2));
            if isempty(ii)
                % reaction has been constrained
                j=floor((I(i,2)-1)/32);
                j=[j+1;I(i,2)-32*j-1];
                j=int8(bitand(R_sub(:,j(1)),uint32(2^j(2)))>0);
                if I(i,3)<0
                    j=-j;
                end;
            else
                % reaction has been left unconstrained
                if i_inc1(I(i,1))==0
                    ii=find(i_A(1,:)==2);
                    t=full(T(I(i,1),i_A(2,ii)));
                    ii(t==0)=[];
                    t(:,t==0)=[];
                    t=t';
                    j=A_sub(:,ii)*t;
                    j((abs(j)./(abs(A_sub(:,ii))*abs(t)))<=param.Tolerance)=0;
                    j=int8(sign(j));
                    i_inc1(I(i,1))=1;
                else
                    j=zeros(m(2)-m(1)+1,1,'int8');
                end;
            end;

            if i_inc(I(i,1))==1
                j1=fm_bin_sub(:,I(i,1));
                i_err1=i_err1|((j.*j1)==-1);
                j1=sign(j1+j);
                fm_bin_sub(:,I(i,1))=j1;
            else
                fm_bin_sub(:,I(i,1))=j;
                i_inc(I(i,1))=1;
            end;
        end;
        clear j j1 R_sub A_sub

        fm_bin_sub(i_err1,:)=-128;

        if iscell(fm_bin)
            fm_bin_sub=fm_bin_sub';
            fwrite_check(fid,fm_bin_sub,'int8');
        else
            fm_bin=fm_bin_sub';
            clear fm_bin_sub
        end;
    end;
else
    if iscell(fm_bin)
        fwrite_check(fid,zeros(0,0,'int8'),'int8');
    end;
end;

if iscell(fm_bin)
    fclose(fid);
    flush_file_cache(fm_bin{1});
end;
if iscell(R)
    fclose(fid1);
    delete(R{1})
end;
if iscell(A)
    fclose(fid2);
    flush_file_cache(A{1});
end;
clearvars -global R;



% save S_unc
buffer_size=get_free_memory;
n_block=min(floor(0.25*buffer_size/(8*A_size(1) + 4*2*(strcmp(pointer_type,'double')+1))),A_size(2));
A_size(1)=A_size(1)-length(i_A_del);

if automatic_splitting|(n_block<A_size(2))|iscell(fm_bin)|iscell(A)
    S_unc={[param.Dir{2},'/S_unc.dat'],A_size,'double'};
    if ~automatic_splitting
        if exist(S_unc{1},'file')
            delete(S_unc{1});
        end;
    end;
else
    S_unc=[];
end;

if iscell(S_unc)
    fid=fopen(S_unc{1},'A');
end;
if iscell(A)
    fid1=fopen(A{1},'r');
end;

if (A_size(1)>0)&iscell(S_unc)
    m=[1-n_block,0];
    while m(2)<R_size(2)
        m=min(m+n_block,R_size(2));

        if iscell(A)
            A_sub=fread(fid1,[A_size(1)+length(i_A_del),m(2)-m(1)+1],'double=>double');
        else
            A_sub=A(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));
        end;
        for i=find(i_A(1,:)==1)
            A_sub(:,i)=A_sub(:,i)/param.CompoundCoefficient(i_A(2,i));
        end;
        A_sub(i_A_del,:)=[];
        if iscell(S_unc)
            fwrite_check(fid,A_sub,'double');
        else
            S_unc=A_sub;
        end;
        clear A_sub
    end;
elseif ~iscell(S_unc)
    A(i_A_del,:)=[];
    S_unc=A;
else
    fwrite_check(fid,[],'double');
end;

if iscell(S_unc)
    fclose(fid);
    flush_file_cache(S_unc{1});
end;
if iscell(A)
    fclose(fid1);
    delete(A{1});
end;
clearvars -global A


i_S_unc=i_A;


if iscell(fm_bin)
    fm_bin=fm_bin{1};
end;
if iscell(S_unc)
    S_unc=S_unc{1};
end;

stats.time=[stats.time,toc];

stats_out=stats;

clearvars -global S_th stats i_A order_A param stage i_perm j switch_j Hamming_weights T_compress



function [n_j,i_current_column,is_zero]=get_next_column

global A stats R order_A i_A param stage j pointer_type




if sum(sum(isinf(order_A),1)==0)==0
    j=zeros(0,1,'int8');
    n_j=[0,0,0];
    i_current_column=[];
    is_zero=1==0;
    return
end;

if iscell(R)
    R_size=R{2};
else
    R_size=size(R);
end;
if R_size(2)>=4294967295
    pointer_type='double';
else
    pointer_type='uint32';
end;

% check if R can be fitted in the memory, and whether there is enough
% memory to contain j and the permutation indices
buffer_size=get_free_memory;
if iscell(R)
    if R_size(2)>param.MaxIntermediates(1)
        j=['Maximum number of intermediate flux modes has been reached (',num2str(R_size(2)),'>',num2str(param.MaxIntermediates(1)),')'];
        n_j=[0,0,0];
        i_current_column=[];
        is_zero=1==0;
        return
    end;
    if (buffer_size-R_size(2)*(4*(strcmp(pointer_type,'double')+1)*2+2+4*R_size(1)))<param.MinAvailableMemory*1048576
        % must be able to store R, 2 permutation arrays and two logical arrays in the memory
        j='Out of memory; cannot fit the current rays and permutation indices in the memory';
        n_j=[0,0,0];
        i_current_column=[];
        is_zero=1==0;
        return
    end
    a=R{1};
    R=matrix_indexing(R);
    buffer_size=buffer_size-4*size(R,1)*size(R,2);
    delete(a)
else
    if R_size(2)>param.MaxIntermediates(1)
        j=['Maximum number of intermediate flux modes has been reached (',num2str(R_size(2)),'>',num2str(param.MaxIntermediates(1)),')'];
        n_j=[0,0,0];
        i_current_column=[];
        is_zero=1==0;
        return
    end;
    if (buffer_size-R_size(2)*(4*(strcmp(pointer_type,'double')+1)*2+2))<param.MinAvailableMemory*1048576
        % must be able to store 2 permutation arrays and two logical arrays in the memory
        j='Out of memory; cannot fit the permutation indices in the memory';
        n_j=[0,0,0];
        i_current_column=[];
        is_zero=1==0;
        return
    end
end;


% calculate number of positive, zero and negative elements in each row of A
n_p0m=zeros(size(i_A,2),3);
if iscell(A)
    A_size=A{2};
else
    A_size=size(A);
end;
mem_block=round(0.25*buffer_size/8);
n_block=floor(mem_block/A_size(1));


if iscell(A)
    fid=fopen(A{1},'r');
end;
m=[1-n_block,0];
while m(2)<A_size(2)
    m=min(m+n_block,A_size(2));
    if iscell(A)
        A_sub=fread(fid,[A_size(1),m(2)-m(1)+1],'double=>double');
    else
        A_sub=A(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));
    end;
    n_p0m(:,1)=n_p0m(:,1)+sum(A_sub>0,2);
    n_p0m(:,2)=n_p0m(:,2)+sum(A_sub==0,2);
    n_p0m(:,3)=n_p0m(:,3)+sum(A_sub<0,2);
    clear A_sub
end
if iscell(A)
    fclose(fid);
    flush_file_cache(A{1});
end;



% criterion for finding the next row
pr=zeros(size(A,1),0);
for i=1:length(param.DynamicOrder)
    if strcmp(param.DynamicOrder{i},'minzero')
        x=n_p0m(:,2);
    elseif strcmp(param.DynamicOrder{i},'maxzero')
        x=-n_p0m(:,2);
    elseif strcmp(param.DynamicOrder{i},'minneg')
        x=n_p0m(:,3);
    elseif strcmp(param.DynamicOrder{i},'maxneg')
        x=-n_p0m(:,3);
    elseif strcmp(param.DynamicOrder{i},'minpos')
        x=n_p0m(:,1);
    elseif strcmp(param.DynamicOrder{i},'maxpos')
        x=-n_p0m(:,1);
    elseif strcmp(param.DynamicOrder{i},'mincomb')
        x=n_p0m(:,1).*n_p0m(:,3);
    elseif strcmp(param.DynamicOrder{i},'maxcomb')
        x=-n_p0m(:,1).*n_p0m(:,3);
    elseif strcmp(param.DynamicOrder{i},'minminposneg')
        x=min(n_p0m(:,[1,3]),[],2);
    else
        error('???')
    end;
    pr=[pr,x];
end;


pr(isinf(order_A(1,:)),:)=Inf;





if sum(~isinf(pr))>0
    [m,i_current_column]=sortrows([order_A',pr]);
    stage=m(1,1);
    i=find((m(:,end)<=(m(1,end)+abs(m(1,end))*0.2))&(sum(m(:,1:(size(m,2)-1))~=m(ones(size(m,1),1),1:(size(m,2)-1)),2)==0));
    i_current_column=i_current_column(i(end));
    n_j=n_p0m(i_current_column,:);
    if i_A(1,i_current_column)==1
        is_zero=1==1;
    else
        is_zero=param.Constraints(i_A(2,i_current_column))==-1;
    end;

    % determine j=int8(sign(A(i_current_column,:)))
    j=zeros(size(A,2),1,'int8');

    mem_block=round(0.25*buffer_size/8);
    n_block=floor(mem_block/A_size(1));


    if iscell(A)
        fid=fopen(A{1},'r');
    end;
    m=[1-n_block,0];
    while m(2)<A_size(2)
        m=min(m+n_block,A_size(2));
        if iscell(A)
            A_sub=fread(fid,[A_size(1),m(2)-m(1)+1],'double=>double');
        else
            A_sub=A(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));
        end;
        j(feval(pointer_type,m(1)):feval(pointer_type,m(2)))=int8(sign(A_sub(i_current_column,:)));
        clear A_sub
    end
    if iscell(A)
        fclose(fid);
        flush_file_cache(A{1});
    end;
else
    j=zeros(0,1,'int8');
    i_current_column=[];
    is_zero=1==0;
    n_j=[0,0,0];
end;



stats.time=[stats.time,toc];
stats.number_FMs=[stats.number_FMs,size(R,2)];
stats.number_constraints=[stats.number_constraints,A_size(1)];
if ~isempty(n_p0m(i_current_column,1).*n_p0m(i_current_column,3))
    stats.number_candidate_FMs=[stats.number_candidate_FMs,[n_p0m(i_current_column,1).*n_p0m(i_current_column,3);0;0]];
else
    stats.number_candidate_FMs=[stats.number_candidate_FMs,[NaN;NaN;NaN]];
end;







function j_comb=find_adjacent_rays(n_j,n_j0)

global S_th i_A order_A R A Hamming_weights param stage stats i_perm j switch_j pointer_type


switch_j=n_j(3)<n_j(1);

if switch_j
    n_j=fliplr(n_j(:)');
end;

buffer_size=get_free_memory;

% internal flag to use harddrive for memory intensive operations and to
% store permutation indices
write_to_harddrive=(buffer_size-size(R,2)*4*(1.5*size(R,1)+5*(strcmp(pointer_type,'double')+1))-param.MinAvailableMemory*1048576)<0;
write_to_harddrive1=(buffer_size-size(R,2)*4*(1.5*size(R,1)+(5+2*8)*(strcmp(pointer_type,'double')+1))-param.MinAvailableMemory*1048576)<0;



th=S_th;
% eliminate rows and columns of S that have not yet been constrained
th(i_A(2,i_A(1,:)==1),:)=[];
i=i_A(2,i_A(1,:)==2);
if (sum(sum(th(:,i)~=0,2)==1)~=length(i))|(sum(sum(th(:,i)~=0,1)==1)~=length(i))
    error('S_th should be in reduced row echelon form . . .')
end;
th(sum(th(:,i)~=0,2)==1,:)=[];
th(:,i_A(2,i_A(1,:)==2))=0;

% eliminate columns that have nonzero elements for the combined rays
R1=0;
m=[1-1048576,0];
nn=length(j);
while m(2)<nn
    m=min(m+1048576,nn);
    jj=feval(pointer_type,find(j(m(1):m(2))~=0)+m(1)-2);
    if ~isempty(jj)
        R1=R1+count_subsets(uint32(ones(size(R,1),1)*(2^32-1)),R,jj,param.MaxThreads);
    end;
end;
R1=R1(1:size(S_th,2),:);
th(:,R1==0)=0;

th(:,sum(th~=0,1)==0)=[];
th=full(th);
th=rank(th)+2;

f_write_error=false;

if exist([param.Dir{1},'/j_comb_cand.dat'],'file')
    delete([param.Dir{1},'/j_comb_cand.dat'])
end;


% display information about the current iteration
if strcmp(param.Display,'iter')
    b=sum(order_A(1,:)==stage);
    a=sum(~isinf(order_A(1,:)));

    if a>1
        str=[sprintf('%3d',a),' constraints left (',sprintf('%2d',b),' in stage ',num2str(stage),')   nr. flux modes: ',sprintf('%5d',size(R,2)),'    nr. candidates: '];
    elseif a==1
        str=[sprintf('%3d',a),' constraint  left (',sprintf('%2d',b),' in stage ',num2str(stage),')   nr. flux modes: ',sprintf('%5d',size(R,2)),'    nr. candidates: '];
    else
        str='';
    end;
    if ~switch_j
        fprintf([str,num2str(n_j(1)*n_j(3)),' (initial: ',num2str(n_j(1)),'/',num2str(n_j(2)),'/',num2str(n_j(3)),'), '])
    else
        fprintf([str,num2str(n_j(1)*n_j(3)),' (initial: ',num2str(n_j(3)),'/',num2str(n_j(2)),'/',num2str(n_j(1)),'), '])
    end;
    str(str~=9)=' ';
end;





% define bit pattern tree parameters
rho=0.5;
n_max=300;
k_max=72;  % should be smaller than the maximum recursion limit in MATLAB of 500

% do not use tree if there are only a few candidate rays
if (n_j(1)*n_j(3)<2*size(R,2))&(min(n_j(1)*n_j(3),size(R,2))<10000)
    n_max=2*size(R,2);
end;



% calculate initial permutation of R
if switch_j
    a=[1,-1,0];
else
    a=[-1,1,0];
end;
n=length(j);
i_perm=ones(1,n,pointer_type);
ni=0;
for a1=a
    m=[1-1048576,0];
    while m(2)<n
        m=min(m+1048576,n);
        jj=feval(pointer_type,find(j(m(1):m(2))==a1)+m(1)-1);
        i_perm((1+ni):(ni+length(jj)))=jj;
        ni=ni+length(jj);
    end;
end;
clear jj

% calculate pattern tree (i_perm - global - is both input and output)
[R_tree,i_leaf,i_node,i_leaf_bound]=get_pattern_tree([1,size(R,2)],rho,n_max,k_max);

% extract + and - tree
i_leaf3=[i_leaf(1,:);i_leaf_bound(1,:)];
i_leaf1=[int64(double(i_leaf_bound(1,:))+1);i_leaf_bound(2,:)];


i_perm1=i_perm;
i_perm=[];

% invert i_perm1
n=length(i_perm1);
i_perm1_inv=ones(1,n,pointer_type);
m=[1-1048576,0];
while m(2)<n
    m=min(m+1048576,n);
    jj=feval(pointer_type,m(1)):feval(pointer_type,m(2));
    i_perm1_inv(i_perm1(jj))=jj;
end;
clear jj



if write_to_harddrive
    R_file=matrix_indexing(R,[param.Dir{1},'/R.dat'],Inf,i_perm1);
    R=[];
    R=matrix_indexing(R_file);
    delete(R_file{1})
    clear R_file
    i_perm1_file=matrix_indexing(i_perm1,[param.Dir{1},'/i_perm1.dat']);
    clear i_perm1
    i_perm1_inv_file=matrix_indexing(i_perm1_inv,[param.Dir{1},'/i_perm1_inv.dat']);
    clear i_perm1_inv
else
    try
        R=R(:,i_perm1);
    catch
        R_file=matrix_indexing(R,[param.Dir{1},'/R.dat'],Inf,i_perm1);
        R=[];
        R=matrix_indexing(R_file);
        delete(R_file{1})
        clear R_file
    end;
end;


j1_inc=zeros(1,n_j(1),'uint32');
n_comb=0;
i_iter=1;

if write_to_harddrive1
    j_comb=[param.Dir{1},'/j_comb.dat'];
    if exist(j_comb,'file')
        delete(j_comb);
    end;
    fid=fopen(j_comb,'W');
    fclose(fid);
else
    j_comb={};
end;


while (sum(j1_inc<size(R_tree,2))>0)&((n_comb+n_j0)<=param.MaxIntermediates(1))&(~f_write_error)
    if (i_iter>1)&strcmp(param.Display,'iter')
        fprintf(', ')
    end;

    % STEP 1: prune candidate rays
    buffer_size=get_free_memory;
    max_blocks=uint32(min(floor(buffer_size*0.35/8388608),512));

    [j_comb_cand,j1_inc,n_comb_cand1]=prune_candidates_rank_pattern_tree_C(R,R_tree,i_node,i_leaf3,i_leaf1,uint16(th),Hamming_weights,max_blocks,j1_inc,[param.Dir{1},'/j_comb_cand.dat'],uint32(param.SSE4),param.MaxThreads,uint32(strcmp(pointer_type,'double')));


    i_comb_cand=NaN;
    if (size(j_comb_cand,1)==1)&(size(j_comb_cand,2)==1)
        flush_file_cache([param.Dir{1},'/j_comb_cand.dat']);
        fid1=fopen([param.Dir{1},'/j_comb_cand.dat'],'r');
        fseek(fid1,0,'eof');
        n_comb_cand=ftell(fid1)/8/(strcmp(pointer_type,'double')+1);
        fseek(fid1,0,'bof');
        try
            j_comb_cand=fread(fid1,[2,n_comb_cand],[pointer_type,'=>',pointer_type]);
        catch
            i_comb_cand=0;
        end;
        fclose(fid1);
        if isnan(i_comb_cand)
            delete([param.Dir{1},'/j_comb_cand.dat']);
        end;
    else
        n_comb_cand=size(j_comb_cand,2);
    end;

    f_write_error=f_write_error|(double(n_comb_cand)~=double(n_comb_cand1));

    if strcmp(param.Display,'iter')
        if (i_iter==1)&(sum(j1_inc==size(R_tree,2))==n_j(1))
            fprintf([num2str(n_comb_cand),' (after pruning), '])
        else
            fprintf([num2str(n_comb_cand),' (after pruning ',sprintf('%3.1f',100*mean(j1_inc)/size(R_tree,2)),'%%), '])
        end;
    end;
    stats.number_candidate_FMs(2,size(stats.number_candidate_FMs,2))=stats.number_candidate_FMs(2,size(stats.number_candidate_FMs,2))+n_comb_cand;


    % STEP 2: check adjacency
    while (i_comb_cand<n_comb_cand)|(isnan(i_comb_cand))&(~f_write_error)
        if isnan(i_comb_cand)
            i_comb_cand=n_comb_cand;
        else
            buffer_size=get_free_memory;
            n_block=min(floor(0.2*buffer_size/8/(strcmp(pointer_type,'double')+1)),n_comb_cand-i_comb_cand);
            fid1=fopen([param.Dir{1},'/j_comb_cand.dat'],'r');
            fseek(fid1,i_comb_cand*8*(strcmp(pointer_type,'double')+1),'bof');
            j_comb_cand=fread(fid1,[2,n_block],[pointer_type,'=>',pointer_type]);
            fclose(fid1);
            if (i_comb_cand>0)&strcmp(param.Display,'iter')
                fprintf(', ')
            end;
            i_comb_cand=i_comb_cand+n_block;
            if n_comb_cand==i_comb_cand
                delete([param.Dir{1},'/j_comb_cand.dat']);
            end;
        end;

        i=adjacency_test_pattern_tree_C(R,j_comb_cand,R_tree,i_node,i_leaf,param.MaxThreads);


        % filter pruned candidate combinations
        i=i~=0;
        try
            j_comb_cand=matrix_indexing(j_comb_cand,'',Inf,i);
        catch
            j_comb_cand=matrix_indexing(j_comb_cand,[param.Dir{1},'/j_comb_cand1.dat'],Inf,i);
            j_comb_cand=matrix_indexing(j_comb_cand);
            delete([param.Dir{1},'/j_comb_cand1.dat'])
        end;
        clear i

        if strcmp(param.Display,'iter')
            fprintf([num2str(size(j_comb_cand,2)),' (after adjacency test)'])
        end;


        n_comb=n_comb+size(j_comb_cand,2);
        if switch_j&(~isempty(j_comb_cand))
            j_comb_cand=j_comb_cand([2,1],:);
        end;
        if write_to_harddrive1
            fid=fopen(j_comb,'A');
            fwrite_check(fid,j_comb_cand,pointer_type);
            fclose(fid);
            flush_file_cache(j_comb);
        else
            try
                j_comb{i_iter}=j_comb_cand;
            catch
                write_to_harddrive1=1==1;
                fid=fopen([param.Dir{1},'/j_comb.dat'],'W');
                for ii=1:length(j_comb)
                    fwrite_check(fid,j_comb{ii},pointer_type);
                end;
                fwrite_check(fid,j_comb_cand,pointer_type);
                j_comb=[param.Dir{1},'/j_comb.dat'];
            end;
        end;
        clear j_comb_cand
        i_iter=i_iter+1;
    end;
end;
if strcmp(param.Display,'iter')
    fprintf(char(10))
end;
if sum(j1_inc<size(R_tree,2))>0
    n_pred=round(n_comb/mean(j1_inc)/size(R_tree,2)+n_j0);
else
    n_pred=0;
end;
clear j1_inc
stats.number_candidate_FMs(3,size(stats.number_candidate_FMs,2))=n_comb;

if write_to_harddrive1
    j_comb={j_comb,[2,n_comb],pointer_type};
else
    j_comb=horzcat(j_comb{:});
end;

if write_to_harddrive
    i_perm1=matrix_indexing(i_perm1_file);
    delete(i_perm1_file{1})
end;

if iscell(j_comb)
    movefile(j_comb{1},[param.Dir{1},'/j_comb1.dat']);
    buffer_size=get_free_memory;
    n_block=floor(0.07*buffer_size/4/(strcmp(pointer_type,'double')+1));

    fid1=fopen([param.Dir{1},'/j_comb1.dat'],'r');
    fid=fopen(j_comb{1},'W');
    fclose(fid);
    x=fread(fid1,n_block,[pointer_type,'=>',pointer_type]);
    while ~isempty(x)
        x=x+1;
        x=i_perm1(x);
        fid=fopen(j_comb{1},'A');
        fwrite_check(fid,x,pointer_type);
        fclose(fid);
        clear x
        flush_file_cache(j_comb{1});
        x=fread(fid1,n_block,[pointer_type,'=>',pointer_type]);
    end;
    fclose(fid1);
    delete([param.Dir{1},'/j_comb1.dat'])
else
    m=[1-1048576,0];
    while m(2)<n_comb
        m=min(m+1048576,n_comb);
        j_comb(:,m(1):m(2))=i_perm1(j_comb(:,m(1):m(2))+1);
    end;
end;

clear i_perm1

if write_to_harddrive
    i_perm1_inv=matrix_indexing(i_perm1_inv_file);
    delete(i_perm1_inv_file{1})
    R_file=matrix_indexing(R,[param.Dir{1},'/R.dat'],Inf,i_perm1_inv);
    R=[];
    clear i_perm1_inv
    R=matrix_indexing(R_file);
    delete(R_file{1})
    clear R_file
else
    try
        R=R(:,i_perm1_inv);
        clear i_perm1_inv
    catch
        R_file=matrix_indexing(R,[param.Dir{1},'/R.dat'],Inf,i_perm1_inv);
        R=[];
        clear i_perm1_inv
        R=matrix_indexing(R_file);
        delete(R_file{1})
        clear R_file
    end;
end;




if n_pred>0
    j_comb=['Maximum number of intermediate flux modes has been reached; predicted number of flux modes after next iteration: ',num2str(n_pred)];
elseif f_write_error
    if exist([param.Dir{1},'/j_comb_cand.dat'],'file')
        delete([param.Dir{1},'/j_comb_cand.dat'])
    end;
    if exist([param.Dir{1},'/j_comb.dat'],'file')
        delete([param.Dir{1},'/j_comb.dat'])
    end;
    stats.number_candidate_FMs([2,3],size(stats.number_candidate_FMs,2))=NaN;
    j_comb='Number of pruned candidate combinations written to the harddrive is incorrect . . .';
end;

if (n_comb+n_j0)>param.MaxIntermediates(1)
    j_comb=['Maximum number of intermediate flux modes has been reached (',num2str(n_comb+n_j0),'>',num2str(param.MaxIntermediates(1)),')'];
end;










function [R_tree_out,i_leaf_out,i_node_out,i_leaf_bound_out]=get_pattern_tree(i_limit,rho_in,n_max_in,k_max_in)

global R rho k n_max k_max S_th i_tree i_leaf i_node i_leaf_bound R_tree_tot i_perm j switch_j

k=0;
rho=rho_in;
n_max=n_max_in;
k_max=k_max_in;
i_tree=2;
R_tree_tot=zeros(size(R,1),ceil((i_limit(2)-i_limit(1)+1)/k_max*1.5),'uint32');
i_leaf=zeros(2,size(R_tree_tot,2));
i_leaf_bound=zeros(2,size(R_tree_tot,2));
i_node=zeros(1,size(R_tree_tot,2));

if (i_limit(2)-i_limit(1)+1)>n_max
    get_pattern_tree_intern(i_limit,zeros(size(S_th,2),1)==1,1);
else
    i_leaf(:,1)=i_limit;
    n_j=full(sparse(1,(2*double(switch_j)-1)*double(j(i_leaf(1,1):i_leaf(2,1)))+2,1,1,3));
    i_leaf_bound(:,1)=[0;n_j(1)]+n_j(3)+i_leaf(1,1)-1;
    i_tree=2;
end;
i_node(1)=i_tree;
i_tree=i_tree-1;

i_leaf=i_leaf(:,1:i_tree);

i_perm1=i_perm;
n=i_limit(1);
for i=1:size(i_leaf,2)
    a=i_leaf(:,i)-i_leaf(1,i)+n;
    i_perm(a(1):a(end))=i_perm1(i_leaf(1,i):i_leaf(2,i));
    i_leaf_bound(:,i)=i_leaf_bound(:,i)-i_leaf(1,i)+n;
    i_leaf(:,i)=a;
    n=n+a(2)-a(1)+1;
end;
clear i_perm1
R_tree_out=R_tree_tot(:,1:i_tree);
i_leaf_out=int64(i_leaf-1);
i_leaf_bound_out=int64(i_leaf_bound(:,1:i_tree)-1);
i_node_out=uint32(i_node(:,1:i_tree)-1);

clear global i_tree
clear global R_tree_tot
clear global i_leaf
clear global i_leaf_bound
clear global i_node
clear global rho
clear global k
clear global n_max
clear global k_max




function get_pattern_tree_intern(i_limit,R_current,i_tree_in)

global R T_compress rho i_perm k n_max k_max i_tree i_leaf i_node i_leaf_bound R_tree_tot j switch_j param pointer_type


k=k+1;

m=[1-1048576,0];
nn=i_limit(2)-i_limit(1)+1;
i_perm1=ones(1,nn,pointer_type);
while m(2)<nn
    m=min(m+1048576,nn);
    i_perm1(m(1):m(2))=i_perm((m(1)+i_limit(1)-1):(m(2)+i_limit(1)-1))-1;
end;

R_tree=R_current;
R_tree_comp=uint32(T_compress*R_tree);

jj=find(~R_tree);

i=[];
a=1;
while (a>rho)
    R_tree(jj(i))=1==1;
    R_tree_comp=uint32(T_compress*R_tree);
    jj(i)=[];
    n=count_subsets(bitcmp(R_tree_comp),R,i_perm1,param.MaxThreads)/(i_limit(2)-i_limit(1)+1);
    n=n(jj);
    a1=a;
    [a,i]=max(n);
end;


if rho-a<a1-rho
    R_tree(jj(i))=1;
end;
R_current1=R_tree;
R_tree=uint32(T_compress*R_tree);

i=is_subset(R_tree,R,i_perm1,param.MaxThreads);
i_perm1=i_perm1+1;

n_i=[];
ni=i_limit(1)-1;
nn=i_limit(2)-i_limit(1)+1;
for aa=[0,1]
    m=[1-1048576,0];
    while m(2)<nn
        m=min(m+1048576,nn);
        jj=find(i(m(1):m(2))==aa)+m(1)-1;
        i_perm((1+ni):(ni+length(jj)))=i_perm1(jj);
        ni=ni+length(jj);
    end;
    n_i=[n_i,ni-i_limit(1)+1];
end;
n_i=[n_i(1),nn-n_i(1)];
clear i i_perm1


if (k<=k_max)&(n_i(2)>n_max)
    i_tree_start=i_tree;
    i_tree=i_tree+1;
    get_pattern_tree_intern([0;(n_i(2)-1)]+i_limit(1)+n_i(1),R_current1,i_tree_start);
    i_node(:,i_tree_start)=i_tree;
else
    R_tree_tot(:,i_tree)=R_tree;
    i_leaf(:,i_tree)=[0;(n_i(2)-1)]+i_limit(1)+n_i(1);
    n_j=full(sparse(1,(2*double(switch_j)-1)*double(j(i_perm(i_leaf(1,i_tree):i_leaf(2,i_tree))))+2,1,1,3));
    i_leaf_bound(:,i_tree)=[0;n_j(1)]+n_j(3)+i_leaf(1,i_tree)-1;
    i_node(:,i_tree)=i_tree+1;
    i_tree=i_tree+1;
end;


if (k<=k_max)&(n_i(1)>n_max)
    get_pattern_tree_intern([0;(n_i(1)-1)]+i_limit(1),R_current,i_tree_in);
else
    R_tree_tot(:,i_tree_in)=uint32(T_compress*R_current);
    i_leaf(:,i_tree_in)=[0;(n_i(1)-1)]+i_limit(1);
    n_j=full(sparse(1,(2*double(switch_j)-1)*double(j(i_perm(i_leaf(1,i_tree_in):i_leaf(2,i_tree_in))))+2,1,1,3));
    i_leaf_bound(:,i_tree_in)=[0;n_j(1)]+n_j(3)+i_leaf(1,i_tree_in)-1;
end;


k=k-1;














function update_DD(n_j,j_comb,i_current_column,is_zero)

global A R i_A order_A param stats j pointer_type


if iscell(j_comb)
    n_comb=j_comb{2}(2);
else
    n_comb=size(j_comb,2);
end;

if is_zero
    n_j_inc=n_j(2);
else
    n_j_inc=n_j(1)+n_j(2);
end;


R_old=R;
buffer_size=get_free_memory;
if (buffer_size<(4*size(R_old,1)+24)*(n_j_inc+n_comb)+param.MinAvailableMemory*1048576)&(~iscell(A))
    % if A is stored in the memory and memory is running low, copy A to the harddrive
    buffer_size=buffer_size+size(A,1)*size(A,2)*8;
    A=matrix_indexing(A,[param.Dir{2},'/A.dat']);
end;
if (buffer_size>4*size(R_old,1)*(n_j_inc+n_comb)+param.MinAvailableMemory*1048576)&(~iscell(A))
    R=ones(size(R_old,1),n_j_inc+n_comb,'uint32');
    R(:)=0;
    buffer_size=buffer_size-4*size(R_old,1)*(n_j_inc+n_comb);
else
    R={[param.Dir{1},'/R.dat'],[size(R_old,1),n_j_inc+n_comb],'uint32'};
end;

R_size=[size(R_old,1),n_j_inc+n_comb];
n_block=floor(0.25*buffer_size/(4*R_size(1)+4*4*(strcmp(pointer_type,'double')+1)));


if iscell(R)
    fid=fopen(R{1},'W');
end;
m=[1-n_block,0];
n=0;
while m(2)<size(R_old,2)
    m=min(m+n_block,size(R_old,2));
    R_sub=R_old(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));

    if is_zero
        j_inc=feval(pointer_type,find(j(feval(pointer_type,m(1)):feval(pointer_type,m(2)))==0));
    else
        j_inc=feval(pointer_type,find(j(feval(pointer_type,m(1)):feval(pointer_type,m(2)))>=0));
    end;
    R_sub=R_sub(:,j_inc);

    if i_A(1,i_current_column)==2
        i=i_A(2,i_current_column)-1;
        i1=floor(i/32);
        i=uint32(2^(i-i1*32));
        i1=i1+1;
        if is_zero
            R_sub(i1,:)=bitand(R_sub(i1,:),bitcmp(i));
        else
            R_sub(i1,:)=bitor(R_sub(i1,:),i*uint32(j(j_inc+m(1)-1,:)'));
        end;
    end;

    if iscell(R)
        fwrite_check(fid,R_sub,'uint32');
    else
        R(:,feval(pointer_type,n+1):feval(pointer_type,n+length(j_inc)))=R_sub;
    end;
    n=n+size(R_sub,2);
    clear R_sub j_inc
end
clear j_inc






if iscell(j_comb)
    fid1=fopen(j_comb{1},'r');
end;

m=[1-n_block,0];
while m(2)<n_comb
    % perform copying in blocks to save memory
    m=min(m+n_block,n_comb);

    if iscell(j_comb)
        j_comb_sub=fread(fid1,[2,m(2)-m(1)+1],[pointer_type,'=>',pointer_type]);
    else
        j_comb_sub=j_comb(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));
    end;
    R_sub=bitor(R_old(:,j_comb_sub(1,:)),R_old(:,j_comb_sub(2,:)));

    if iscell(R)
        fwrite_check(fid,R_sub,'uint32');
    else
        R(:,feval(pointer_type,n_j_inc+m(1)):feval(pointer_type,n_j_inc+m(2)))=R_sub;
    end;
    clear R_sub j_comb_sub
end;
clear R_old

if iscell(R)
    fclose(fid);
    flush_file_cache(R{1});
end;
if iscell(j_comb)
    fclose(fid1);
    flush_file_cache(j_comb{1});
end;


stats.id=[stats.id,i_A(:,i_current_column)];


i_new=[1:size(i_A,2)]~=i_current_column;
i_A=i_A(:,i_new);
order_A=order_A(:,i_new);

A_old=A;
if iscell(A_old)
    A_size=[A_old{2}(1)-1,n_j_inc+n_comb];
else
    A_size=[size(A_old,1)-1,n_j_inc+n_comb];
end;


if ((buffer_size>8*A_size(1)*A_size(2)+param.MinAvailableMemory*1048576)&(~iscell(A_old))&(~iscell(R)))|(A_size(1)==0)
    A=ones(A_size(1),A_size(2),'double');
    A(:)=0;
    buffer_size=buffer_size-8*A_size(1)*A_size(2);
else
    A={[param.Dir{2},'/A.dat'],A_size,'double'};
    if iscell(A_old)
        if disk_space(param.Dir{1})>(prod(A_old{2})*8+1e10)
            % copy A_old to fast medium if enough space is available
            movefile(A_old{1},[param.Dir{1},'/A_old.dat'])
            A_old{1}=[param.Dir{1},'/A_old.dat'];
        else
            movefile(A_old{1},[param.Dir{2},'/A_old.dat'])
            A_old{1}=[param.Dir{2},'/A_old.dat'];
        end;
    end;
end;
if iscell(A_old)
    n_block=floor(0.25*buffer_size/(4*8*(A_size(1)+1)+4*6*(strcmp(pointer_type,'double')+1)));
else
    n_block=floor(0.25*buffer_size/(3*8*(A_size(1)+1)+4*4*(strcmp(pointer_type,'double')+1)));
end;

tol_range=[0,10.^[-16:0],Inf];
t=([tol_range]'<0)+0;
tol1=[];
a_max=0;

if A_size(1)>0
    if iscell(A)
        fid=fopen(A{1},'W');
    end;
    if iscell(A_old)
        fid1=fopen(A_old{1},'r');
        n1=A_old{2}(2);
    else
        n1=size(A_old,2);
    end;
    m=[1-n_block,0];
    n=0;
    while m(2)<n1
        m=min(m+n_block,n1);
        if iscell(A_old)
            A_sub=fread(fid1,[A_size(1)+1,m(2)-m(1)+1],'double=>double');
        else
            A_sub=A_old(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));
        end;

        if is_zero
            j_inc=feval(pointer_type,find(j(feval(pointer_type,m(1)):feval(pointer_type,m(2)))==0));
        else
            j_inc=feval(pointer_type,find(j(feval(pointer_type,m(1)):feval(pointer_type,m(2)))>=0));
        end;
        A_sub=A_sub(i_new,j_inc);

        if iscell(A)
            fwrite_check(fid,A_sub,'double');
        else
            A(:,feval(pointer_type,n+1):feval(pointer_type,n+length(j_inc)))=A_sub;
        end;
        n=n+size(A_sub,2);
        clear A_sub j_inc
    end

    % empty global variable j
    j=[];

    if iscell(j_comb)
        fid2=fopen(j_comb{1},'r');
    end;

    m=[1-n_block,0];
    while m(2)<n_comb
        % perform copying in blocks to save memory
        m=min(m+n_block,n_comb);

        if iscell(j_comb)
            j_comb_sub=fread(fid2,[2,m(2)-m(1)+1],[pointer_type,'=>',pointer_type]);
        else
            j_comb_sub=j_comb(:,feval(pointer_type,m(1)):feval(pointer_type,m(2)));
        end;
        if iscell(A_old)
            fseek(fid1,0,'bof');
            j_comb_sub1=j_comb_sub-1;
            j_comb_sub1=mod(j_comb_sub1,2*n_block);
            j_comb_sub=j_comb_sub-j_comb_sub1;
            j_comb_sub1=j_comb_sub1+1;
            A1_sub=ones(A_size(1)+1,m(2)-m(1)+1,'double');
            A2_sub=ones(A_size(1)+1,m(2)-m(1)+1,'double');
            A1_sub(:)=0;
            A2_sub(:)=0;
            m1=[1-2*n_block,0];
            while m1(2)<n1
                m1=min(m1+2*n_block,n1);
                X=fread(fid1,[A_size(1)+1,m1(2)-m1(1)+1],'double=>double');
                jj=j_comb_sub(1,:)==feval(pointer_type,m1(1));
                A1_sub(:,jj)=X(:,j_comb_sub1(1,jj));
                jj=j_comb_sub(2,:)==feval(pointer_type,m1(1));
                A2_sub(:,jj)=X(:,j_comb_sub1(2,jj));
                clear X jj
            end;
            A1_sub=A1_sub';
            A2_sub=A2_sub';
        else
            A1_sub=A_old(:,j_comb_sub(1,:))';
            A2_sub=A_old(:,j_comb_sub(2,:))';
            j_comb_sub1=[];
        end;
        clear j_comb_sub j_comb_sub1

        c1=-A2_sub(:,i_current_column);
        c2=A1_sub(:,i_current_column);
        A1_sub=A1_sub(:,i_new);
        A2_sub=A2_sub(:,i_new);
        if param.all_integer
            g=gcd_fast(c1,c2);
            i=(g>1)&(~isnan(g));
            c1(i,:)=round(c1(i,:)./g(i,:));
            c2(i,:)=round(c2(i,:)./g(i,:));
            clear g i
        end;
        A1_sub=A1_sub.*c1(:,ones(1,size(A1_sub,2)));
        A2_sub=A2_sub.*c2(:,ones(1,size(A2_sub,2)));
        clear c1 c2
        B=abs(A1_sub)+abs(A2_sub);
        A1_sub=A1_sub+A2_sub;
        A2_sub=B;
        clear B

        if param.all_integer
            % condition rows in A1_sub with the gcd
            A2_sub=abs(A1_sub);
            A2_sub(A2_sub==0)=NaN;
            a_max=max(a_max,max(max(A2_sub)));
            i=min(A2_sub,[],2)>=10;
            g=1./gcd_columns(A1_sub(i,:)')';
            A1_sub(i,:)=round(A1_sub(i,:).*g(:,ones(1,size(A1_sub,2))));
        else
            % assign 0 to elements for which the relative value <= tol
            A2_sub=abs(A1_sub)./A2_sub;
            t=t+histc(A2_sub(:)',tol_range)';
            A1_sub(A2_sub<=param.Tolerance)=0;
            A2_sub=A2_sub((A2_sub>param.Tolerance)&(A2_sub<=100*param.Tolerance));
            if ~isempty(A2_sub)
                tol1=[tol1;min(min(A2_sub));max(max(A2_sub))];
            end;

            % condition rows in A1_sub
            A2_sub=abs(A1_sub);
            A2_sub(A2_sub==0)=NaN;
            a_max=max(a_max,max(max(A2_sub)));
            g=min(A2_sub,[],2);
            i=(~isnan(g))&((g<1e-3)|(g>1e3));
            g=2.^round(-log(g)/log(2));
            A1_sub(i,:)=A1_sub(i,:).*g(i,ones(1,size(A1_sub,2)));
        end;

        A1_sub=A1_sub';
        if iscell(A)
            fwrite_check(fid,A1_sub,'double');
        else
            A(:,feval(pointer_type,n_j_inc+m(1)):feval(pointer_type,n_j_inc+m(2)))=A1_sub;
        end;
        clear A1_sub A2_sub j_comb_sub
    end;


    if iscell(A)
        fclose(fid);
        flush_file_cache(A{1})
    end;
    if iscell(A_old)
        fclose(fid1);
    end;
    if iscell(j_comb)
        fclose(fid2);
    end;
end;


if iscell(j_comb)
    delete(j_comb{1})
end;
if iscell(A_old)
    delete(A_old{1})
end;
clear A_old j_comb

if ~isempty(tol1)
    a=['Warning: combination of rays gave relative coefficients in the range ',num2str(min(tol1)),' ',num2str(max(tol1)),' . . .'];
    disp(a)
    stats.warnings=[stats.warnings;{a}];
end;

t(end)=[];
if sum(t)>0
    t=t/sum(t);
end;

stats.relative_coefficient_distribution=[stats.relative_coefficient_distribution,t];
stats.max_coefficient=[stats.max_coefficient,a_max];
if (a_max>2^53)&param.all_integer
    if ~strcmp(param.Display,'off')
        disp('Integer precision lost because coefficients larger than 2^53; switching to floating point precision . . .')
    end;
    param.all_integer=1==0;
end;



function initialise_DD(A_in,i_A_in,order_A_in,R_in,S_th_in,param_in,i_A_resolved)

global A S_th R T_compress stats i_A order_A Hamming_weights param stage i_perm

A=[];
R=[];
i_perm=[];

i=uint16([0:65535]');
Hamming_weights=zeros(65536,1,'uint8');
for j=1:16
    Hamming_weights=Hamming_weights+uint8(bitget(i,j));
end;

if size(A_in,1)~=size(R_in,2)
    error('Number of rows in A should be equal to number of columns in R . . .')
end;


stats=struct('max_coefficient',ones(1,size(i_A_resolved,2))*max([max(A_in(:)'),0]),'time',ones(1,size(i_A_resolved,2))*toc,'number_FMs',ones(1,size(i_A_resolved,2))*size(R_in,2),'number_constraints',ones(1,size(i_A_resolved,2))*size(A_in,2),'number_candidate_FMs',zeros(3,size(i_A_resolved,2)),'id',i_A_resolved,'relative_coefficient_distribution',zeros(18,size(i_A_resolved,2)),'warnings',{cell(0,1)});


i_A=i_A_in;
S_th=S_th_in;
order_A=order_A_in;
param=param_in;
stage=min(order_A(1,:));

A=A_in';
i=[1:size(S_th,2)];
T_compress=full(sparse(floor([i-1]/32)+1,i,2.^mod(i-1,32)));
if (mod(size(T_compress,1),2)==1)&(size(T_compress,1)>=5)
    % enforce the MEX files to use 64 bit operations
    T_compress=[T_compress;zeros(1,size(T_compress,2))];
end;
R=uint32(T_compress*double(R_in));









function [C_blocked,rev_change,C_compression,p_compression]=LP_compression(S_in,rev,i_excl,tol,LPsolver)


if nargin<5
    LPsolver='';
end;
if nargin<4
    tol=1e-10;
end;
if nargin<3
    i_excl=[];
end;

if ~islogical(i_excl)
    i_excl=full(sparse(1,i_excl,1,1,size(S_in,2)))>0;
else
    i_excl=i_excl(:)';
end;

rev_change=zeros(1,size(S_in,2));

i_S=find((sum(S_in~=0,1)>0)&(rev~=-1));
S=S_in(:,i_S);
S=condition_matrix(S,'Round2power','on');


[a,i]=sortrows([-rev(:,i_S)',sum(S~=0,1)',sum(abs(S),1)']);
if sum(sum(S~=round(S)))==0
    [S_rref,p,i]=rref_fast(S,'ColumnOrder',i,'Integer','on');
else
    [S_rref,p,i]=rref_fast(S,'ColumnOrder',i);
end;
S_rref(sum(S_rref~=0,2)==0,:)=[];
S=S(i(1:size(S_rref,1)),:);

S_rev=S_rref(1:sum(p&(rev(:,i_S)==1)),:)';
S_irr=S_rref(sum(p&(rev(:,i_S)==1))+(1:sum(p&(rev(:,i_S)==0))),:)';

% detect irreversible cycles
i=sum(S_rev((p==0)&(rev(:,i_S)==1),:)~=0,1)>0;
rev_cycle=(sum(S_rev(:,i)~=0,2)>0)'&(rev(:,i_S)==1);

i_blocked=[];
j=find(rev(:,i_S)==0);

% find blocked irreversible reactions
for i=j
    [x,f_opt,exitflag]=LP_solve(zeros(size(S_irr,2),1),-S_irr(j(j~=i),:),zeros(sum(j~=i),1),S_irr(i,:),1,'Integer','off','Solver',LPsolver);
    if exitflag==1
        i_blocked=[i_blocked,i_S(i)];
    end;
end;

% find constrained/blocked reversible reactions
k=find((rev(:,i_S)==1)&(~rev_cycle));
for i=k
    ii=find(S_rev(i,:)~=0);
    if length(ii)==1
        s_rev=S_rev(:,ii)*sign(S_rev(i,ii));
    else
        error('???')
    end;
    [x1,f_opt1,exitflag1]=LP_solve(zeros(size(S_irr,2),1),-S_irr(j,:),-s_rev(j,:),[],[],'Integer','off','Solver',LPsolver);
    [x2,f_opt2,exitflag2]=LP_solve(zeros(size(S_irr,2),1),-S_irr(j,:),s_rev(j,:),[],[],'Integer','off','Solver',LPsolver);


    if (exitflag1==1)&(exitflag2==1)
        i_blocked=[i_blocked,i_S(i)];
    elseif exitflag1==1
        rev_change(i_S(i))=1;
    elseif exitflag2==1
        rev_change(i_S(i))=-1;
    end;
end;

C_blocked=full(sparse(1:length(i_blocked),i_blocked,1,length(i_blocked),size(S_in,2)));



if nargout>2
    % perform full LP compression

    % add constraints and re-perform Gaussian elimination
    S=S_in(:,i_S);
    clear S_rref S_rev S_irr
    S(:,sum(C_blocked(:,i_S),1)>0)=[];
    i_S(sum(C_blocked(:,i_S),1)>0)=[];
    S(:,rev_change(i_S)==-1)=-S(:,rev_change(i_S)==-1);
    rev(rev_change~=0)=0;
    [a,i]=sortrows([-rev(:,i_S)',sum(S~=0,1)',sum(abs(S),1)']);
    if sum(sum(S~=round(S)))==0
        [S_rref,p,i]=rref_fast(S,'ColumnOrder',i,'Integer','on');
    else
        [S_rref,p,i]=rref_fast(S,'ColumnOrder',i,'Integer','off');
    end;
    S_rref(sum(S_rref~=0,2)==0,:)=[];

    %all_integer=full(sum(sum(S_rref~=round(S_rref))))==0;

    S_rev=S_rref(1:sum(p&(rev(:,i_S)==1)),:)';
    S_irr=S_rref(sum(p&(rev(:,i_S)==1))+(1:sum(p&(rev(:,i_S)==0))),:);

    %s=ones(1,size(S_irr,2));
    %if all_integer
    %    [S_irr,x,s]=condition_matrix_integer(S_irr);
    %end;
    S_irr=S_irr';

    [S_irr1,s1,s2]=condition_matrix(S_irr,'Round2power','on');
    s1=s1(:);
    s2=s2(:);

    i=sum(S_rev((p==0)&(rev(:,i_S)==1),:)~=0,1)>0;
    rev_cycle=(sum(S_rev(:,i)~=0,2)>0)'&(rev(:,i_S)==1);

    C=zeros(0,length(i_S));
    i_zero=zeros(1,length(i_S))==1;

    % find compression for irreversible reactions
    j=find(rev(:,i_S)==0);
    for i=j(:,~i_excl(i_S(j)))
        i_eq=[j(:,(j~=i)&(i_zero(:,j))),i];
        i_ineq=j((j~=i)&(~i_zero(:,j)));


        [x,f_opt,exitflag,output]=LP_solve(-S_irr1(j(j~=i),:)'*(1./(s1(j(j~=i),:))),-S_irr1(i_ineq,:),zeros(length(i_ineq),1),S_irr1(i_eq,:),[zeros(length(i_eq)-1,1);-1],'Integer','off','Solver',LPsolver,'ConditionCoeffs',{s2(:),s1(i_ineq,:),s1(i_eq,:)});
        if exitflag==1
            x=find(output.b~=0);
            if length(x)~=1
                error('???')
            end;
            a=1./abs(output.b(x));
            if a>1
                output.A(x,:)=a*output.A(x,:);
                output.b(x,:)=a*output.b(x,:);
            end;
            %if all_integer&((sum(sum(output.A~=round(output.A)))>0)|(sum(sum(output.b~=round(output.b)))>0))
            %    error('Could not find integer solution A b . . .')
            %end;
            [X,Y]=rref_mult(output.A,output.b,S_irr,[],tol);
            if iscell(Y)
                x=full([Y{1}';Y{2}']);  %full([Y{1}'.*s;Y{2}']);
                g=gcd_columns(x);
                x(:,~isnan(g))=x(:,~isnan(g))./g([1;1],~isnan(g));
                if sum(x(:,2)~=round(x(:,2)))==0
                    l=lcm_columns(x(2,:)');
                else
                    l=1e100;
                end;
                if l<2^53
                    x=x(1,:).*round(l./x(2,:));
                else
                    x=x(1,:)./x(2,:);
                end;
            else
                x=Y'; %.*s;
            end;
            x(i_zero)=0;

            [a,l]=round_fraction(x(:,x~=0)');
            if l<1e6
                x=round(l*x);
            end;
            g=gcd_columns(x(:,x~=0)');
            if ~isnan(g)
                x=x./g;
            end;

            ii=find(sum(C(:,i)~=0)>0);
            if ~isempty(ii)
                y=C(ii,:)*x(i)-C(ii,i)*x;
                y(:,i)=0;
                g=gcd_columns(y')';
                y(~isnan(g),:)=y(~isnan(g),:)./g(~isnan(g),ones(1,size(y,2)));
                C(ii,:)=y;
            end;
            C=[C;x];
            i_zero(i)=1==1;
        end;
    end;


    % find compression for reversible reactions
    R=condition_matrix(S_rev((rev(:,i_S)~=1)|(~p),:),'Round2power','on');
    R=((1./mean(abs(R),2))*ones(1,size(R,2))).*R;
    R(isnan(R))=0;
    R=abs(eye(size(R,2))+1-(sparse(diag(1./diag(R'*R)))*((R'*R).^2))./(ones(size(R,2),size(R,1))*(R.^2)))<tol;
    k=find(rev(:,i_S)==1);
    for i=k(:,(~i_excl(:,i_S(k)))&p(:,k))
        ii=find(S_rev(i,:)~=0);
        if length(ii)==1
            s_rev=S_rev(:,ii)*sign(S_rev(i,ii));
        else
            error('???')
        end;

        x=[];
        if ((sum(s_rev(k)~=0)==2)&(sum(s_rev(j)~=0)==0)) | ((sum(s_rev(j)<0)==0)&(sum(s_rev(k)~=0)==1)) | ((sum(s_rev(j)>0)==0)&(sum(s_rev(k)~=0)==1))
            x=s_rev';

        elseif sum(R(:,ii))>0
            ii1=find(R(:,ii));
            ii1=ii1(1);
            if sum(sum(S_rev(:,[ii,ii1])~=round(S_rev(:,[ii,ii1]))))==0 %all_integer
                x=rref_fast(S_rev(:,[ii,ii1])','ColumnOrder',find(rev(:,i_S)==0),'Integer','on');
            else
                x=rref_fast(S_rev(:,[ii,ii1])','ColumnOrder',find(rev(:,i_S)==0),'Integer','off');
            end;
            x=x(2,:);

        end;

        if ~isempty(x)
            x(i_zero)=0;
            [a,l]=round_fraction(x(:,x~=0)');
            if (l<1e6)
                x=round(l*x);
            end;
            g=gcd_columns(x(:,x~=0)');
            if ~isnan(g)
                x=x./g;
            end;

            for i1=1:size(C,1)
                i2=find((C(i1,:)~=0)&i_zero);
                if x(i2)~=0
                    x=x*C(i1,i2)-C(i1,:)*x(i2);
                    g=gcd_columns(x(:,x~=0)');
                    if ~isnan(g)
                        x=x./g;
                    else
                        x=x/C(i1,i2);
                    end;
                end;
            end;

            ii=find(sum(C(:,i)~=0)>0);
            if ~isempty(ii)
                y=C(ii,:)*x(i)-C(ii,i)*x;
                y(:,i)=0;
                g=gcd_columns(y')';
                y(~isnan(g),:)=y(~isnan(g),:)./g(~isnan(g),ones(1,size(y,2)));
                C(ii,:)=y;
            end;

            C=[C;x];
            i_zero(i)=1==1;
            R((i_zero(:,k)*S_rev(k,:))>0,:)=0==1;
        end;
    end;


    C_compression=zeros(size(C,1),size(S_in,2));
    C_compression(:,i_S)=C;
    p_compression=zeros(1,size(S_in,2));
    p_compression(:,i_S)=i_zero;
    C_compression(:,rev_change==-1)=-C_compression(:,rev_change==-1);
    clear c i_zero

    if sum(sum(C_compression~=round(C_compression)))==0
        [C_compression,p_compression]=rref_fast(C_compression,'ColumnOrder',find(p_compression),'Integer','on');
    else
        [C_compression,p_compression]=rref_fast(C_compression,'ColumnOrder',find(p_compression));
    end;
    g=gcd_columns(C_compression')';
    C_compression(~isnan(g),:)=C_compression(~isnan(g),:)./g(~isnan(g),ones(1,size(C_compression,2)));

    % check potential inaccuracies in solution in case of floating points
    a=abs(C_compression);
    a(a==0)=NaN;
    i=find((max(a,[],2)./min(a,[],2))>1e6);
    if (~isempty(i))&(sum(sum(C_compression~=round(C_compression)))>0)
        p_compression((sum(C_compression(i,:)~=0,1)>0)&(p_compression==1))=0==1;
        C_compression(i,:)=[];
    end;
end;





function Y=matrix_indexing(X,Y_filename,i1,i2,varargin)


param=function_arguments(struct('FlushFileCache','on','BufferSize',5e7),varargin);

if strcmp(lower(param.FlushFileCache),'on')
    param.FlushFileCache=1==1;
elseif ~strcmp(lower(param.FlushFileCache),'off')
    error(['Unknown option for FlushFileCache: ',param.FlushFileCache])
else
    param.FlushFileCache=1==0;
end;

strclass={'double','single','int8','int16','int32','int64','uint8','uint16','uint32','uint64'};
byteclass=[8,4,1,2,4,8,1,2,4,8];


if nargin<4
    i2=Inf;
end;
if nargin<3
    i1=Inf;
end;
if nargin<2
    Y_filename='';
end;

if iscell(X)
    if length(X)~=3
        error('X has invalid format - should be {''filename'', matrix size, ''data type''} ');
    end;
    if ~ischar(X{1})
        error('X has invalid format - should be {''filename'', matrix size, ''data type''} ');
    end;
    if exist(X{1},'file')==0
        error('filename X does not exist ');
    end;
    X_filename=X{1};
    n=dir(X{1});
    n=n.bytes;

    if ~ischar(X{3})
        error('X has invalid format - should be {''filename'', matrix size, ''data type''} ');
    end;
    b=find(strcmp(strclass,X{3}));
    if length(b)~=1
        error(['Invalid data class for X: ',X{3}])
    end;
    b=byteclass(b);
    X_type=X{3};

    if ~isnumeric(X{2})
        error('X has invalid format - should be {''filename'', matrix size, ''data type''} ');
    elseif length(X{2})>2
        error('X has invalid format - should be {''filename'', matrix size, ''data type''} ');
    end;
    X_size=X{2};
    if length(X_size)==1
        X_size=[X_size,1];
    end;

    if prod(X_size(:))*b~=n
        error(['Matrix size (',num2str(X_size(1)),'*',num2str(X_size(2)),'*',num2str(b),') and file size (',num2str(n),') are not equal'])
    end;

    X_mem=1==0;

else
    X_type=class(X);
    b=find(strcmp(strclass,X_type));
    if length(b)~=1
        error(['Invalid data class for X: ',X_type])
    end;
    X_size=size(X);
    X_mem=1==1;
end;

if isempty(Y_filename)
    Y_mem=1==1;
elseif ischar(Y_filename)
    Y_mem=1==0;
else
    error('Second argument should be empty (Y stored in memory) or a filename (Y stored on harddrive) ')
end;

% determine if all indices are taken for 1st and 2nd dimension
i1_all=0==1;
if length(i1)==1
    i1_all=isinf(i1);
end;
i2_all=0==1;
if length(i2)==1
    i2_all=isinf(i2);
end;

Y_size=[0,0];
if i1_all
    Y_size(1)=X_size(1);
else
    if islogical(i1)
        if length(i1)~=X_size(1)
            error('size of i1 and X are not equal')
        end;
        Y_size(1)=sum(i1);
    elseif isnumeric(i1)
        if (min(i1)<1)|(max(i1)>X_size(1))
            error(['values of i1 must lie within the range [1..',num2str(X_size(1)),']'])
        end;
        Y_size(1)=length(i1);
    else
        error('i1 should be logical or numeric index')
    end;
end;
if i2_all
    Y_size(2)=X_size(2);
else
    if islogical(i2)
        if length(i2)~=X_size(2)
            error('size of i2 and X are not equal')
        end;
        Y_size(2)=sum(i2);
    elseif isnumeric(i2)
        if (min(i2)<1)|(max(i2)>X_size(2))
            error(['values of i2 must lie within the range [1..',num2str(X_size(2)),']'])
        end;
        Y_size(2)=length(i2);
    else
        error('i2 should be logical or numeric index')
    end;
end;

Y_type=X_type;

if Y_mem
    Y=ones(Y_size(1),Y_size(2),Y_type);
    Y(:)=0;
else
    Y={Y_filename,Y_size,Y_type};
end;

if (~X_mem)&(~Y_mem)
    if strcmp(X_filename,Y_filename)
        error('Filenames must be different . . .')
    end
end;

if (X_mem)&(Y_mem)
    mem_block=round(param.BufferSize/byteclass(strcmp(strclass,X_type)));
    if i1_all
        n_block=floor(mem_block/X_size(1));
    else
        n_block=floor(mem_block/length(i1));
    end;
    if i2_all
        n=X_size(2);
    else
        n=length(i2);
    end;

    if n_block>=n
        if (i1_all)&(i2_all)
            Y(:)=X+0;
        elseif (i1_all)&(~i2_all)
            Y=X(:,i2);
        elseif (~i1_all)&(i2_all)
            Y=X(i1,:);
        elseif (~i1_all)&(~i2_all)
            Y=X(i1,i2);
        end;

    elseif n_block==0

        if i1_all
            n1=X_size(1);
        else
            n1=length(i1);
        end;
        nY=[0,0];
        for i=1:n
            m=[1-mem_block,0];
            if ~i2_all
                if islogical(i2)
                    nY(2)=nY(2)+double(i2(i));
                end;
            end;
            nY(1)=0;
            while m(2)<n1
                m=min(m+mem_block,n1);

                j=m(1):m(2);
                j1=j;
                if ~i1_all
                    if islogical(i1)
                        j=find(i1(j))+m(1)-1;
                        j1=(1:length(j))+nY(1);
                        nY(1)=nY(1)+length(j);
                    else
                        j=i1(j);
                    end;
                end;

                if ~isempty(j)
                    if i2_all
                        Y(j1,i)=X(j,i);
                    elseif ~islogical(i2)
                        Y(j1,i)=X(j,i2(i));
                    elseif i2(i)
                        Y(j1,nY(2))=X(j,i);
                    end;
                end;
            end;
        end;

    else
        m=[1-n_block,0];
        nY=[0,0];
        while m(2)<n
            m=min(m+n_block,n);

            j=m(1):m(2);
            j2=j;
            if ~i2_all
                if islogical(i2)
                    j=find(i2(j))+m(1)-1;
                    j2=(1:length(j))+nY(2);
                    nY(2)=nY(2)+length(j);
                else
                    j=i2(j);
                end;
            end;

            if ~isempty(j)
                if i1_all
                    Y(:,j2)=X(:,j);
                else
                    Y(:,j2)=X(i1,j);
                end;
            end;
        end;
    end;

elseif (~X_mem)&(Y_mem)
    mem_block=round(param.BufferSize/byteclass(strcmp(strclass,X_type)));
    n_block=floor(mem_block/X_size(1));
    n=X_size(2);

    if n_block>=n
        clear Y;
        fid=fopen(X_filename,'r');
        y=fread(fid,X_size,[X_type,'=>',X_type]);
        fclose(fid);
        if (i1_all)&(i2_all)
            Y=y;
        elseif (i1_all)&(~i2_all)
            Y=y(:,i2);
        elseif (~i1_all)&(i2_all)
            Y=y(i1,:);
        elseif (~i1_all)&(~i2_all)
            Y=y(i1,i2);
        end;

    elseif n_block==0
        fid=fopen(X_filename,'r');
        nY=[0,0];
        if i1_all
            n1=X_size(1);
        else
            n1=length(i1);
            if ~islogical(i1)
                i1=sortrows([[1:length(i1)]',i1],2);
            end;
        end;
        for i=1:n
            m=[1-mem_block,0];
            nY(1)=0;
            if ~i2_all
                if islogical(i2)
                    nY(2)=nY(2)+double(i2(i));
                    no_read=~i2(i);
                else
                    if sum(i2==i)==0
                        nY(2)=NaN;
                    else
                        nY(2)=find(i2==i);
                    end;
                    no_read=sum(i2==i)==0;
                end;
            else
                no_read=1==0;
            end;
            if ~no_read
                while m(2)<n1
                    m=min(m+mem_block,n1);
                    y=fread(fid,m(2)-m(1)+1,[X_type,'=>',X_type]);

                    j=m(1):m(2);
                    if ~i1_all
                        if islogical(i1)
                            j1=(find(i1(j)));
                            j=(1+nY(1)):(length(j1)+nY(1));
                            nY(1)=nY(1)+sum(i1((m(1)):(m(2))));
                        else
                            j=(i1(:,2)>=m(1))&(i1(:,2)<=m(2));
                            j1=i1(j,2)-m(1)+1;
                            j=i1(j,1);
                        end;
                    else
                        j1=j-m(1)+1;
                    end;



                    if ~isempty(j)
                        if i2_all
                            Y(j,i)=y(j1);
                        elseif ~islogical(i2)
                            if ~isnan(nY(2))
                                Y(j,nY(2))=y(j1);
                            end;
                        elseif i2(i)
                            Y(j,nY(2))=y(j1);
                        end;
                    end;
                end;
            else
                fseek(fid,X_size(1)*byteclass(strcmp(strclass,X_type)),'cof');
            end;
        end;
        fclose(fid);
    else
        fid=fopen(X_filename,'r');
        m=[1-n_block,0];
        nY=[0,0];
        if ~i2_all
            if ~islogical(i2)
                i2=sortrows([[(1):(length(i2))]',(i2)],2);
            end;
        end;
        while m(2)<n
            m=min(m+n_block,n);
            y=fread(fid,[X_size(1),m(2)-m(1)+1],[X_type,'=>',X_type]);

            j=(m(1)):(m(2));
            if ~i2_all
                if islogical(i2)
                    j1=(find(i2(j)));
                    j=(1+nY(2)):(length(j1)+nY(2));
                    nY(2)=nY(2)+sum(i2((m(1)):(m(2))));
                else
                    j=(i2(:,2)>=m(1))&(i2(:,2)<=m(2));
                    j1=i2(j,2)-m(1)+1;
                    j=i2(j,1);
                end;
            else
                j1=j-m(1)+1;
            end;

            if ~isempty(j)
                if i1_all
                    Y(:,j)=y(:,j1);
                else
                    Y(:,j)=y(i1,j1);
                end;
            end;
        end;
        fclose(fid);
    end;

    if param.FlushFileCache
        flush_file_cache(X_filename);
    end;

elseif (X_mem)&(~Y_mem)
    if exist(Y_filename,'file')
        delete(Y_filename)
    end;

    mem_block=round(param.BufferSize/byteclass(strcmp(strclass,X_type)));
    if i1_all
        n_block=floor(mem_block/X_size(1));
    else
        n_block=floor(mem_block/length(i1));
    end;
    if i2_all
        n=X_size(2);
    else
        n=length(i2);
    end;

    if n_block>=n
        fid=fopen(Y_filename,'W');
        if (i1_all)&(i2_all)
            fwrite_check(fid,X,X_type);
        elseif (i1_all)&(~i2_all)
            fwrite_check(fid,X(:,i2),X_type);
        elseif (~i1_all)&(i2_all)
            fwrite_check(fid,X(i1,:),X_type);
        elseif (~i1_all)&(~i2_all)
            fwrite_check(fid,X(i1,i2),X_type);
        end;
        fclose(fid);

    elseif n_block==0
        fid=fopen(Y_filename,'W');
        if i1_all
            n1=X_size(1);
        else
            n1=length(i1);
        end;
        for i=1:n
            m=[1-mem_block,0];
            while m(2)<n1
                m=min(m+mem_block,n1);

                j=(m(1)):(m(2));
                if ~i1_all
                    if islogical(i1)
                        j=(find(i1(j))+m(1)-1);
                    else
                        j=i1(j);
                    end;
                end;

                if ~isempty(j)
                    if i2_all
                        fwrite_check(fid,X(j,i),X_type);
                    elseif ~islogical(i2)
                        fwrite_check(fid,X(j,i2(i)),X_type);
                    elseif i2(i)
                        fwrite_check(fid,X(j,i),X_type);
                    end;
                end;
            end;
        end;
        fclose(fid);

    else
        fid=fopen(Y_filename,'W');
        m=[1-n_block,0];
        while m(2)<n
            m=min(m+n_block,n);

            j=(m(1)):(m(2));
            if ~i2_all
                if islogical(i2)
                    j=(find(i2(j))+m(1)-1);
                else
                    j=i2(j);
                end;
            end;

            if ~isempty(j)
                if i1_all
                    fwrite_check(fid,X(:,j),X_type);
                else
                    fwrite_check(fid,X(i1,j),X_type);
                end;
            end;
        end;
        fclose(fid);
    end;

    if param.FlushFileCache
        flush_file_cache(Y_filename);
    end;

else
    if ((~i1_all)&(~islogical(i1)))|((~i2_all)&(~islogical(i2)))
        error('Numerical indices are not supported for harddrive-harddrive indexing . . .')
    end;
    mem_block=round(param.BufferSize/byteclass(strcmp(strclass,X_type)));
    n_block=floor(mem_block/X_size(1));
    n=X_size(2);

    if n_block>=n
        fid=fopen(X_filename,'r');
        y=fread(fid,X_size,[X_type,'=>',X_type]);
        fclose(fid);

        fid1=fopen(Y_filename,'W');
        if (i1_all)&(i2_all)
            fwrite_check(fid1,y,Y_type);
        elseif (i1_all)&(~i2_all)
            fwrite_check(fid1,y(:,i2),Y_type);
        elseif (~i1_all)&(i2_all)
            fwrite_check(fid1,y(i1,:),Y_type);
        elseif (~i1_all)&(~i2_all)
            fwrite_check(fid1,y(i1,i2),Y_type);
        end;
        fclose(fid1);

    elseif n_block==0
        fid=fopen(X_filename,'r');
        fid1=fopen(Y_filename,'W');
        if i1_all
            n1=X_size(1);
        else
            n1=length(i1);
            if ~islogical(i1)
                i1=sortrows([[(1):(length(i1))]',(i1)],2);
            end;
        end;
        for i=1:n
            m=[1-mem_block,0];
            if ~i2_all
                if islogical(i2)
                    no_read=~i2(i);
                else
                    error('???')
                end;
            else
                no_read=1==0;
            end;
            if ~no_read
                while m(2)<n1
                    m=min(m+mem_block,n1);
                    y=fread(fid,m(2)-m(1)+1,[X_type,'=>',X_type]);

                    if ~i1_all
                        if islogical(i1)
                            j1=(find(i1((m(1)):(m(2)))));
                        else
                            error('???')
                        end;
                    else
                        j1=(1):(m(2)-m(1)+1);
                    end;



                    if ~isempty(j1)
                        if i2_all
                            fwrite_check(fid1,y(j1),Y_type);
                        elseif ~islogical(i2)
                            error('???')
                        elseif i2(i)
                            fwrite_check(fid1,y(j1),Y_type);
                        end;
                    end;
                end;
            else
                fseek(fid,X_size(1)*byteclass(strcmp(strclass,X_type)),'cof');
            end;
        end;
        fclose(fid1);
        fclose(fid);

    else
        fid=fopen(X_filename,'r');
        fid1=fopen(Y_filename,'W');
        m=[1-n_block,0];

        if ~i2_all
            if ~islogical(i2)
                error('???');
            end;
        end;
        while m(2)<n
            m=min(m+n_block,n);
            y=fread(fid,[X_size(1),m(2)-m(1)+1],[X_type,'=>',X_type]);

            if ~i2_all
                if islogical(i2)
                    j1=(find(i2((m(1)):(m(2)))));
                else
                    error('???');
                end;
            else
                j1=(1):(m(2)-m(1)+1);
            end;

            if ~isempty(j1)
                if i1_all
                    fwrite_check(fid1,y(:,j1),Y_type);
                else
                    fwrite_check(fid1,y(i1,j1),Y_type);
                end;
            end;
        end;
        fclose(fid1);
        fclose(fid);
    end;

    if param.FlushFileCache
        flush_file_cache(X_filename);
        flush_file_cache(Y_filename);
    end;
end;





function n=get_free_memory

global param

x=memory_OS_indep;
n=x.MaxPossibleArrayBytes;

% determine amount of memory that is used by all MATLAB variables and the
% remaining amount of memory, in case an upper limit is defined
a=whos;
a=sum(vertcat(a.bytes));
b=whos('global');
a=a+sum(vertcat(b.bytes));
n=min(n,max(0,param.MaxAvailableMemory*1048576-a));


function fwrite_check(fid,x,mode)

n=fwrite(fid,x,mode);

if n~=prod(size(x))
    filename=fopen(fid);
    x=x(:);
    n1=length(x);
    x(1:n)=[];
    n2=n;
    while ~isempty(x)
        a=input(['Warning: could only write ',num2str(n2),' of ',num2str(n1),' elements to ',filename,'; press "Enter" to retry or "q" to quit . . .'],'s');
        if strcmp(lower(a),'q')
            fclose(fid);
            error(['Could not write ',filename,' . . .'])
        end;
        n=fwrite(fid,x,mode);
        x(1:n)=[];
        n2=n2+n;
    end;
end;
