function param=function_arguments(param_in,param_in_raw)
% Convert input arguments from cell (varargin) to struct. param_in gives
% the expected input parameters with default values, param_in_raw should be
% passed varargin.
%
%  Copyright (c) 2015, Jan Bert van Klinken
%

param=param_in;
if round(length(param_in_raw)/2)*2~=length(param_in_raw)
    error('Invalid syntax of option parameters . . .')
end;
for i=1:round(length(param_in_raw)/2)
    if ischar(param_in_raw{2*i-1})
        if isfield(param,param_in_raw{2*i-1})
            param.(param_in_raw{2*i-1})=param_in_raw{2*i};
        else
            error(['Unknown option: ',param_in_raw{2*i-1},' . . .']);
        end;
    else
        disp(param_in_raw{2*i-1})
        error(['Illegal parameter assignment . . .'])
    end;
end;
