function [st,si]=sms_33

%% Returns a matrix representing the set of equalities that define
%% a 3x3 semi-magic square, i.e. constant row and column sums.

st=[1 1 1 -1 -1 -1 0  0  0;
1 1 1 0  0  0  -1 -1 -1;
0 1 1 -1  0  0  -1  0  0;
1 0 1 0  -1  0  0  -1  0;
1 1 0 0  0  -1  0  0  -1];

si=ones(1,9);
