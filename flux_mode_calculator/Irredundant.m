% This function remove redundancy, i.e. supersets, in a given monotone Boolean function and returns an irredundant monotone Boolean function.
% MBF is a monotone Boolean function as a Binary matrix. The type of MBF, i.e. CNF or DNF, does not matter during the procedure of removing redundancy.

function MBF = Irredundant(MBF)

nrow = size(MBF, 1);
chk = false(1, nrow); % chk(i) is one if MBF(i,:) is a superset and should be removed

% Sort rows in ascending order based on the number of elements in each clause/monomial
summation = sum(MBF, 2);
[~,I] = sort(summation);
MBF = MBF(I, :);

for i=1:(nrow-1)
    a = MBF(i, :);
    for j=(i+1):nrow
        if ~chk(j) % if clause/monomial j is not already checked for being removed
            chk(j) = isequal(and(MBF(j,:), a), a);
%             chk(j) = ~ any(and(MBF(j,:), a)==a);
%             temp = rowfun(@eqchk, table(MBF((i+1):nrow,:), repmat(a, nrow-i, 1)));
%             chk((i+1):nrow) = table2array(temp)' & chk((i+1):nrow);
        end
    end
end
MBF(chk,:) = [];


% i = 1;
% while( ~isempty(MBF) && i <= size(MBF,1))
%     a = MBF(i, :);
%     chk = false(1, size(MBF,1));
%     for j=(i+1):size(MBF,1)
%             chk(j) = isequal(and(MBF(j,:), a), a);
%     end
%     MBF(chk,:) = [];
%     i = i+1;
% end


return
end


