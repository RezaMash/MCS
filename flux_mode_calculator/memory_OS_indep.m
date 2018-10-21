function [m1,m2]=memory_OS_indep
% OS independent version of MATLAB's 'memory' command; unix part needs to
% be extended

if ispc
    % if OS is windows
    [m1,m2]=memory;

elseif isunix
    % if OS is unix
    [a,b]=unix('free  | grep Mem');
    stats=str2double(regexp(b,'[0-9]*','match'));
    % uncertain whether MATLAB has access to the full amount of free memory; multiply by 0.9 to be sure
    m1.MaxPossibleArrayBytes=round(sum(stats([3,5,6]))*0.9)*1024;
    m1.MemAvailableAllArrays=round(sum(stats([3,5,6]))*0.9)*1024;
    m1.MemUsedMATLAB=NaN;
    
    m2.PhysicalMemory=struct('Total',stats(1)*1024,'Available',sum(stats([3,5,6]))*1024);
    
else
    error(['Unknown operating system: ',computer])
end;