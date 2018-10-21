addpath('../CellNetAnalyzer/');
[S, rev, transposed_null] = read(model_name, reduced);
[m,n] = size(S);
I = eye(n);
T = zeros([1, n]);
T(target_reaction) = 1;

Iirrev = zeros(n:0)
for i = 1:n
   if rev(i) == 0   % for irreversible reactions only
       Iirrev = [Iirrev, (I(:,i))];
   end
end

D = [I, -T', S'];
rev_dual = [ones([1, n]),[0], ones([1,m])]

tic;
[efm_bin,S_unc,id_unc]=calculate_flux_modes(D,rev_dual);
fileID = fopen(strcat('logs/',log_file), 'a+');
fprintf(fileID, '%d,%10.3f,', target_reaction, toc);

D_EFM = efm_bin';
fileID = fopen('to_send_to_java.txt', 'w');
for i = 1:size(efm_bin,2)

    if D_EFM(i, n+1) ~=0
      support = abs(sign(D_EFM(i,1:n)));
      positive_support = (abs(sign(D_EFM(i,1:n))) + sign(D_EFM(i,1:n)))*0.5;
      I_support = (support .* int8(rev)) + positive_support .* int8(ones([1, n]) - rev);
      if isequal(I_support,T) || isequal(I_support,zeros([1,n]))
        continue;
      end
      fprintf(fileID, '%d', I_support);
      fprintf(fileID, '\n');
    end
end
exit
