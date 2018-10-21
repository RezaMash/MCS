function r = index_nzeros_rev_adj(A, rev) %indices of non-zero reversible and positive of irreversible elements
r = [];
for i = 1:size(A,2)
    if ((rev(i) ~= 0 && A(i) ~=0) || (A(i) > 0))
        r = [r,i];
    end
end
