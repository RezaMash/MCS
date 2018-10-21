%%

clc
clear


% generate random stoichiometric matrix
n=[10 20];
S=[diag(floor(rand(n(1),1)*3)+1),round(0.7*randn(n(1),n(2)-n(1)))];
A=1;
while rank(A)<n(1)
    A=round(randn(n(1)*1.5,n(1))*0.3);
end;
S=A*S;
rev=rand(1,n(2))>0.6;
% S = [-1,1,-1,-1,0;0,0,1,1,-1];
% rev = [1,0,0,0,0];
%--------------------------------------------------------------
sz = size(S);
m = sz(1);
n = sz(2);
I = eye(n);
D = [S', I];
for i = 1:n
   if rev(i) == 0   % for irreversible reactions only
       D = [D, -(I(:,i))];
   end
end
t = input('Target reaction(s) index:');
T = zeros([1 n]);
for i = 1:n
    if ismember(i, t)
        T(i) = -1;
    end
end
D = [D,T']

sz = size(D);
new_m = sz(1)
new_n = sz(2)



[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(D,ones([1, new_n]));

D_EFM = efm_bin'

% % % % show help of both functions
% help calculate_flux_modes
% help bin2num_flux_modes
% minimal_size = 10000;
% for i = 1:n
%     nzeros = index_nzeros(D_EFM(i, m+1:m+n));
%     w = D(i, end);
%     if w ~= 0 && minimal_size > size(nzeros,2)
%         minimal_size = size(nzeros,2)
%     end
% end
sz = size(D_EFM);

for i = 1:sz(1)
    nzeros = index_nzeros(D_EFM(i, m+1:m+n));
    w = D_EFM(i, end);
    if w ~= 0 % && minimal_size == size(nzeros,2)
        disp('Found MCS:')
        disp(nzeros)
        disp('corresponded EM for dual:')
        disp(D_EFM(i,:))
    end
end

% X = [1,1,0,0,2;0,1,2,3,4];
% disp(index_nzeros(X(1,2:5)));

% calculate binary EFMs
% [efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(S,rev);
% disp('EFM:')
% disp(efm_bin)
% disp('S_unc:')
% disp(S_unc)
% disp('id_unc:')
% disp(id_unc)
% disp('T:')
% disp(T)
% disp('stats:')
% disp(stats)

% efm_bin;
% calculate coefficients of binary EFMs and check for consistency
% [efm,err]=bin2num_flux_modes(efm_bin,S);
