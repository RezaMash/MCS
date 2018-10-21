%
% diary('the_result.txt')
[S, rev, transposed_null] = read(model_name, reduced);

% rev = [0,0
sz = size(S);
m = sz(1);
n = sz(2);
rev = zeros([1,n]);
I = eye(n);
num_irrev = 0;
t = input('Target reaction(s) index:');
T = zeros([1 n]);
for i = 1:n
    if ismember(i, t)
        T(i) = 1;
    end
end

Iirrev = zeros(n:0)
num_irrev = 0;
for i = 1:n
   if rev(i) == 0   % for irreversible reactions only
       Iirrev = [Iirrev, (I(:,i))];
       num_irrev = num_irrev + 1;
   end
end

% D = [I, -T', S'];
D = [I, -T', -Iirrev, S'];
disp(full(D));

sz = size(D);
new_m = sz(1)
new_n = sz(2)

% rev = [ones([1, n]),[0], ones([1,m])]
rev = [ones([1, n]),zeros([1, n+1]), ones([1,m])]

[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(D,rev);

D_EFM = efm_bin';

sz = size(D_EFM);
number_of_MCSs = 0;

for i = 1:size(efm_bin,2)
    nzeros = index_nzeros(D_EFM(i,1:n));
    w = D_EFM(i, n+1);
    if w ~= 0
        number_of_MCSs = number_of_MCSs + 1;
        disp(['Found MCS ', num2str(number_of_MCSs),':']);
        disp(nzeros);
        disp('corresponded EM for dual:')
        disp(D_EFM(i,:))
    end
end

number_of_MCSs
