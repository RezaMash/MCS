function r = index_negatives(A) %index of non zero elements in an array
r = [];
for i = 1:size(A,2)
    if A(i) < 0
        r = [r, i];
    end
end
