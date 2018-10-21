addpath('./flux_mode_calculator/');
[S, rev, transposed_null] = read(model_name, reduced);

[rows, columns] = size(S);
[rows_tn, columns_tn] = size(transposed_null);
tic;
[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(transposed_null, ones([1,columns_tn]));
efm_bin = efm_bin';


if 0 == exist('logs', 'dir')
    mkdir('logs');
end

if target_reaction == -1
  if 0 == exist('temp', 'dir')
      mkdir('temp');
  end
  for target = 1:size(efm_bin,2)
    candidate_MCSs = efm_bin(efm_bin(:,target)>0,:);
    fileID = fopen(strcat('./temp/to_send_to_java_', num2str(target), '.txt'), 'w');
    for candidate_MCS = candidate_MCSs.'
      support = abs(sign(candidate_MCS'));
      negative_support = (abs(sign(candidate_MCS')) - sign(candidate_MCS'))*0.5;
      I_coordianted_support = (support .* int8(rev)) + negative_support .* int8(ones([1, columns]) - rev);
      fprintf(fileID, '%d', I_coordianted_support);
      fprintf(fileID, '\n');
    end
  end
else
  candidate_MCSs = efm_bin(efm_bin(:,target_reaction)>0,:);
  fileID = fopen('to_send_to_java.txt', 'w');
  for candidate_MCS = candidate_MCSs.'
    support = abs(sign(candidate_MCS'));
    negative_support = (abs(sign(candidate_MCS')) - sign(candidate_MCS'))*0.5;
    I_coordianted_support = (support .* int8(rev)) + negative_support .* int8(ones([1, columns]) - rev);
    fprintf(fileID, '%d', I_coordianted_support);
    fprintf(fileID, '\n');
  end
end
fileID = fopen(strcat('logs/',log_file), 'a+');
fprintf(fileID, 'EFMs of transposed null, %10.3f , -\n', toc);
