
if reduced
  S = dlmread(strcat('./models/', model_name, '/reduced_matrix_floated_dlm_matlab_format.txt'))
  rev = dlmread(strcat('./models/', model_name, '/reduced_reversibility_dlm_matlab_format.txt'))
  transposed_null = dlmread(strcat('./models/', model_name, '/reduced_null_basis_floated_dlm_matlab_format.txt'))
else
  S = dlmread(strcat('./models/', model_name, '/full_matrix_floated_dlm_matlab_format.txt'))
  rev = dlmread(strcat('./models/', model_name, '/reversibility_dlm_matlab_format.txt'))
  transposed_null = dlmread(strcat('./models/', model_name, '/null_basis_floated_dlm_matlab_format.txt'))
end

sz = size(S);
rows = sz(1);
columns = sz(2);


% transposed_null = null(S)';
sz = size(transposed_null);
rows_tn = sz(1);
columns_tn = sz(2);

all_null_rev = ones([1,columns_tn])

addpath('./flux_mode_calculator/');
[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(transposed_null, all_null_rev);

D_EFM = efm_bin';
number_of_MCSs = 0;

% if 0 == exist('logs', 'dir')
%     mkdir('logs');
% end

% if 0 == exist(strcat('./logs/',model_name))
%     mkdir('logs', model_name);
% end
%
% if reduced
%     fileID = fopen(strcat('./logs/',model_name,'_reduced_matrix_and_minimal_supports_in_row_space.txt'), 'w');
% else
%     fileID = fopen(strcat('./logs/',model_name,'_matrix_and_minimal_supports_in_row_space.txt'), 'w');
% end
% fprintf(fileID, 'The reversibility status:\n');
% fprintf(fileID, repmat('%d\t', 1, columns), rev);
% if reduced
%     fprintf(fileID, '\nThe reduced stoichiometric matrix:\n');
% else
%     fprintf(fileID, '\nThe stoichiometric matrix:\n');
% end
% for i = 1:rows
%     fprintf(fileID, repmat('%d\t', 1, columns), S(i,:));
%     fprintf(fileID, '\n');
% end
%
% fprintf(fileID, '\n\nThe vectors with minimal support in row-space:\n');
%
% for i = 1:size(efm_bin,2)
%       efm = D_EFM(i,1:columns_tn);
%       fprintf(fileID, repmat('%d\t', 1, columns), efm);
%       fprintf(fileID, '\n');
% end
% display ('finished printing efms\n')
if reduced
    fileID = fopen(strcat('./logs/',model_name,'_mcs_statics_reduced.txt'), 'w');
else
    fileID = fopen(strcat('./logs/',model_name,'_mcs_statics.txt'), 'w');
end
header = 'T \tAll \tReds \tper \tSupers \tper \tTotal percentage of wastes';
formatSpec = '\n%d \t%d \t%d \t%4.2f \t%d \t%4.2f \t%4.2f';
Data = zeros(0, 7);
fprintf(fileID, strcat(model_name, '\nreactions(irreversibles): %d(%d)\nmetabolites: %d\n'), columns, columns- nnz(rev), rows);
fprintf(fileID, header);
for target_reaction = 1: columns
    MCSs = zeros(0,columns);
    all_results = 0;
    for i = 1:size(efm_bin,2)
        efm = D_EFM(i,1:columns_tn);
        if efm(1,target_reaction) > 0
            efm;
            support = abs(sign(efm));
            negative_support = (abs(sign(efm)) - sign(efm))*0.5;
            I_coordianted_support = (support .* int8(rev)) + negative_support .* int8(ones([1, columns]) - rev);
            MCSs = [MCSs; I_coordianted_support];
            all_results = all_results + 1;
            % index_nzeros (I_coordianted_support)
        end
    end

    sortrows(MCSs);
    reduntants = 0;
    supersets = 0;
    for i = 1:size(MCSs)
        for j = 1:(i-1)
            if isequal(MCSs(i, :), MCSs(j, :))
                reduntants = reduntants + 1;
                break
            end
            if isequal(MCSs(i, :).*MCSs(j, :), MCSs(i, :))
                supersets = supersets + 1;
                break
            end
            if isequal(MCSs(i, :).*MCSs(j, :), MCSs(j, :))
                supersets = supersets + 1;
                break
            end
        end
    end
    Data = [Data; 10000-10000*(supersets+reduntants)/all_results, target_reaction, all_results, reduntants, 10000*reduntants/all_results, supersets, 10000*supersets/all_results];
end
Data = sortrows(Data);
for i = 1: columns
    fprintf(fileID, formatSpec, Data(i,2), Data(i,3), Data(i,4), double(Data(i,5))/100, Data(i,6), double(Data(i,7))/100, (10000-double(Data(i,1)))/100);
end
