function [small_stoic small_signs]=small_init

small_signs=[1   1   1   0   0   1   1   0   1   1   0];

small_stoic=[ 1  0  0 -1  0 -1  0  1  0  1  0;
0  0  0  0 -2  2  0 -2  0  0  0;
4  0  0  2  2  4 -3  0  0  0  0;
0 -1  0  0  0  0  1  0  0  0  0;
-1  0  0  0  0  0  0  0  1  0  0;
0  0  1  0  0  0  0  0  0 -1  0;
0  0  0  -0.5  -0.5  0  0  0  0  0  1];
