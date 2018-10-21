addpath('./FK/src');
addpath('./flux_mode_calculator/');
addpath('../CellNetAnalyzer/');
[S, rev, transposed_null] = read(model_name, reduced);
tic;
[ems,S_unc,id_unc,T,stats]=calculate_flux_modes(S, rev);
fileID = fopen(strcat('logs/',log_file), 'a+');
fprintf(fileID, 'EFMs Finding, %10.3f , -\n', toc);
tic;
ems= ems';
if target_reaction == -1
  for target = 1:size(ems,2)
    ems_temp = ems;
    ems_temp(ems_temp ~= 0) = 1;
    ems_temp = ems_temp(ems_temp(:,target)~=0,:);
    ems_temp(:,target) = 0;
    mat = Preprocessing( ems_temp );
    mat = Irredundant(mat);

    [IMFK_cnf, IMFK_cpu, Imcnf_len]= Improved_MFK_Dualization( mat );
    fprintf(fileID, 'MCSs for target: %d,%10.3f,%d\n',target, toc, size(IMFK_cnf,1));
  end
else
  tic;
  ems(ems ~= 0) = 1 ;
  ems = ems(ems(:,target_reaction)~=0,:);
  ems(:,target_reaction) = 0;

  mat = ems;
  mat = Preprocessing( mat );
  mat = Irredundant(mat);

  [IMFK_cnf, IMFK_cpu, Imcnf_len]= Improved_MFK_Dualization( mat );
  fprintf(fileID, 'MCSs for target: %d,%10.3f,%d\n',target_reaction, toc, size(IMFK_cnf,1));

end
exit

exit
