addpath('./flux_mode_calculator/'); 
addpath('../CellNetAnalyzer/')
startcna
startcna(1);
[S, rev, transposed_null] = read(model_name, reduced);
tic;
[ems,S_unc,id_unc,T,stats]=calculate_flux_modes(S, rev);
ems= ems';
fileID = fopen(strcat('logs/',log_file), 'w');
fprintf(fileID, 'Time for the main process: %10.3f seconds\n', toc);
tic;
ems = ems(ems(:,target_reaction)~=0,:);
ems(:,target_reaction) = 0;
MCSs = CNAcomputeCutsets(ems);
fprintf(fileID, 'Time for the dualization process: %10.3f seconds\n', toc);
fprintf(fileID, 'Number of MCSs: %d', size(MCSs,1))

exit
