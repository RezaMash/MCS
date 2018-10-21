S = dlmread(strcat('./models/', model_name, '/full_matrix_floated_dlm_matlab_format.txt'))
rev = dlmread(strcat('./models/', model_name, '/reversibility_dlm_matlab_format.txt'))
[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(S,rev);
