%%


% diary('the_result.txt')



% sbml = sbmlimport('/home/reza/Workspace/MetabolicNetworks/Data_for_Dualization/Data/BIOMD0000000034.xml');
% sbml.Compartments;
[S, rev, transposed_null] = read(model_name, reduced);

% S = M
%
% S = [1,-1,0,0,0,0,0,0,0,0,0,0;0,0,1,-1,0,-1,1,0,0,0,0,0;0,-1,0,1,-1,0,0,0,0,0,0,0;0,1,0,0,1,0,0,1,0,-1,0,0;0,0,0,0,0,1,-1,-1,-1,0,0,0;0,0,0,0,0,0,0,0,1,0,-1,0;0,0,0,0,0,0,0,0,0,1,0,-1]
%
% disp(full(S));
% rev = [1,0,0,0,0];

% rev = [0,0
[m,n] = size(S);
I = eye(n);
% num_irrev = 0;
% t = input('Target reaction(s) index:');

T = zeros([1, n]);
T(target_reaction) = 1;

Iirrev = zeros(n:0)
% num_irrev = 0;
for i = 1:n
   if rev(i) == 0   % for irreversible reactions only
       Iirrev = [Iirrev, (I(:,i))];
       % num_irrev = num_irrev + 1;
   end
end

D = [I, -T', S'];
% D = [I, -T', -Iirrev, S'];


rev_dual = [ones([1, n]),[0], ones([1,m])]
% rev = [ones([1, n]),zeros([1, n+1]), ones([1,m])]
tic;
[efm_bin,S_unc,id_unc]=calculate_flux_modes(D,rev_dual);
fileID = fopen(strcat('logs/',log_file), 'w');
fprintf(fileID, 'Time for the main process: %10.3f seconds\n', toc);

D_EFM = efm_bin';
fileID = fopen('to_send_to_java.txt', 'w');
for i = 1:size(efm_bin,2)

    if D_EFM(i, n+1) ~=0
      support = abs(sign(D_EFM(i,1:n)));
      positive_support = (abs(sign(D_EFM(i,1:n))) + sign(D_EFM(i,1:n)))*0.5;
      I_support = (support .* int8(rev)) + positive_support .* int8(ones([1, n]) - rev);
    % positives = index_nzeros_rev_adj(D_EFM(i,1:n),rev);

    % w = D_EFM(i, n+1);
    % if w ~= 0
        % to_print = [];
        % if rev(target_reaction) ~= 0
        %   to_print = I_support;
        % else
        %   to_print = support;
        % end
        if isequal(I_support,T) || isequal(I_support,zeros([1,n]))
          continue;
        end
        fprintf(fileID, '%d', I_support);
        fprintf(fileID, '\n');
        % number_of_MCSs = number_of_MCSs + 1;
        % disp(['Found MCS ', num2str(number_of_MCSs),':']);
        % disp(positives);
        % disp('corresponded EM for dual:');
        % disp(D_EFM(i,:));
    end
end
% disp(rev);
% number_of_MCSs
exit
