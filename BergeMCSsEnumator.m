addpath('./flux_mode_calculator/');
addpath('../CellNetAnalyzer/')
startcna
startcna(1);
[S, rev, transposed_null] = read(model_name, reduced);
tic;
[ems,S_unc,id_unc,T,stats]=calculate_flux_modes(S, rev);
ems= ems';
fileID = fopen(strcat('logs/',log_file), 'a+');
fprintf(fileID, 'EFMs Finding, %10.3f , -\n', toc);
if target_reaction == -1
  for target = 1:size(ems,2)
    tic;
    ems_temp = ems(ems(:,target)~=0,:);
    ems_temp(:,target) = 0;
    MCSs = CNAcomputeCutsets(ems_temp);
    fprintf(fileID, 'MCSs for target: %d,%10.3f,%d\n',target, toc, size(MCSs,1));
  end
else
  tic;
  ems = ems(ems(:,target_reaction)~=0,:);
  ems(:,target_reaction) = 0;
  MCSs = CNAcomputeCutsets(ems);
  fprintf(fileID, 'MCSs for target: %d,%10.3f,%d\n',target_reaction, toc, size(MCSs,1));
end
exit
