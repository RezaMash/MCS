%%


% diary('the_result.txt')
null

if null_is ready == false
    if sbml == true
        sbml = sbmlimport(strcat('./models/', model_name, '/sbml.xml'));
        [S,objSpecies,objReactions] = sbml.getstoichmatrix;
        S = full(S);
    else
        S = dlmread(strcat('./models/', model_name, '/matlab_dlm.txt'))
    end
else
end
% sbml.Compartments;


% S = [1,-1,0,0,0,0,0,0,0,0,0,0;0,0,1,-1,0,-1,1,0,0,0,0,0;0,-1,0,1,-1,0,0,0,0,0,0,0;0,1,0,0,1,0,0,1,0,-1,0,0;0,0,0,0,0,1,-1,-1,-1,0,0,0;0,0,0,0,0,0,0,0,1,0,-1,0;0,0,0,0,0,0,0,0,0,1,0,-1]

N = null(S)'
sz = size(N);
m = sz(1);
n = sz(2);

rev = ones([1,n])
t = [1];
addpath('./flux_mode_calculator/')
[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(N,rev);


D_EFM = efm_bin';

sz = size(D_EFM);
number_of_MCSs = 0;

if 0 == exist('logs', 'dir')
    mkdir('logs')
end

if 0 == exist(strcat('./logs/',model_name))
    mkdir('logs', model_name)
end
fileID = fopen(strcat('./logs/',model_name,'/mcs_of_',string(t),'.txt'), 'w')
MCSs = zeros(0,n);

for i = 1:size(efm_bin,2)
    % D_EFM(i,1:n)
    % rev

    % nzeros = index_nzeros(D_EFM(i,1:n))
    efm = D_EFM(i,1:n);
    if efm(1,t) > 0
        % negatives = index_negatives(efm);
        % number_of_MCSs = number_of_MCSs + 1;
        mcs = (abs(sign(efm)) - sign(efm))*0.5;
        MCSs = [MCSs;mcs];
        % disp(['Found MCS ', num2str(number_of_MCSs),':']);
        % disp(negatives);
        % disp(efm);
    end
end

sortrows(MCSs);
for i = 1:size(MCSs)
    display(index_nzeros(MCSs(i,1:n)))
    fprintf(fileID,'%d ', index_nzeros(MCSs(i,1:n)));
    fprintf(fileID,'\n');
end
% number_of_MCSs

% exit()
