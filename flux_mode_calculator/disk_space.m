function n=disk_space(d)
% get amount of disk space in bytes


n=0;
if ispc
    [x,y]=dos(['dir ',d]);
    y=split_str(y,char(10));
    while isempty(y{end})
        y(end)=[];
    end;    
    y=lower(strtrim(y{end}));
    y=split_str(y,' ');
    i=find(strcmp(y,'bytes'));
    if length(i)==1
        if strcmp(y{i+1},'free')
            y=y{i-1};
            y(y=='.')='';
            n=str2double(y);
            if isnan(n)
                n=0;
            end;
        end;
    end;
    
elseif isunix
    [err,n]=unix(['df --block-size=1 ',d,' | tail -1 | awk ''{print $4}'' ']);
    n=str2double(n);
    if err~=0
        n=0;
    end;
    
else
    error('Unknown operating system . . .')
end;