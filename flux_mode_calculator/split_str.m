function str_out=split_str(str,delim)
% split string 'str' for occurrences of char 'delim'

if iscell(str)
    str=char(str);
end;

str_out=cell(size(str,1),max(sum(str==delim,2))+1);
for i=1:size(str,1)
    x=str(i,:);
    j=[0,find(x==delim),length(x)+1];
    for i1=2:length(j)
        str_out{i,i1-1}=strtrim(x((j(i1-1)+1):(j(i1)-1)));
    end;
end;