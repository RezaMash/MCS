function [st, si] = sms_55

%% Returns a matrix "st" representing the set of equalities that define
%% a 5x5 semi-magic square, i.e. constant row and column sums.
%% Also returns, "si" a vector of ones (nonnegativity for all variables).
%%
%% This is a test case for joint generation of a hypergraph & its dual.

st=[1 1 1 1 1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 1 1 1 1 0 0 0 0 0 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0;
1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1 0 0 0 0 0;
1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1;
0 1 1 1 1 -1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0;
1 0 1 1 1 0 -1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1 0 0 0;
1 1 0 1 1 0 0 -1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1 0 0;
1 1 1 0 1 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1 0;
1 1 1 1 0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1];

si=ones(1,25);
