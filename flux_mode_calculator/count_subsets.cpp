/*

This mex function is part of the MATLAB function CALCULATE_FLUX_MODES

N = count_subsets(R_subset,R,i_perm,max_threads)

with as input 
    R_subset                - bit pattern 
    R                       - the binary (non-zero) flux modes
    i_perm                  - indices of R over which to count the bit patterns
    max_threads             - maximum number of threads to use

and as output
    N                       - N[i] is the number of flux modes in R(:,i_perm) that have x as subset,
                              with x the complement of R_subset and with the i-th bit set to one 


  Copyright (c) 2015, Jan Bert van Klinken

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include "mex.h"
#include "stdio.h"
#include "math.h"
#include "omp.h"
#include "stdint.h"

#define VAR_INIT0(RSIZE,ISIZE)                                                        \
        const uint32_t m=(uint32_t)mxGetM(prhs[0])*32/RSIZE;                    \
        uint##RSIZE##_t* const R = (uint##RSIZE##_t*) mxGetPr(prhs[1]);         \
        uint##RSIZE##_t* const R_subset = (uint##RSIZE##_t*) mxGetPr(prhs[0]);  \
        ISIZE* const i_perm = (ISIZE*) mxGetPr(prhs[2]);

#define VAR_INIT(RSIZE,ISIZE)                                                   \
            uint64_t j1;                                                        \
            int64_t i1;                                                         \
            uint32_t i2,i3,j2;                                                  \
            uint##ISIZE##_t *N_local = (uint##ISIZE##_t*) malloc(m*64*sizeof(uint##ISIZE##_t));                 \
            uint##RSIZE##_t R1;                                                                                 \
            for (i2=0;i2<m*64;i2++) {N_local[i2]=0;}
            

# define FOR_LOOP_32(ISIZE)                                                                             \
            for (i1=0;i1<ny;i1++) {                                                                           \
                j1=(uint64_t) i_perm[i1];                                                             \
                j1*=m;                                                                                        \
                i3=0;                                                                                         \
                while ((( R_subset[i3]| R[j1+i3])==0xffffffff)&&(i3<m)) {i3++;}                       \
                if (i3==m) {                                                                                  \
                    R1=R[j1];                                        \
                    N_local[0]+=(uint##ISIZE##_t)(( 0x00000001 & R1)!=0);  \
                    N_local[1]+=(uint##ISIZE##_t)(( 0x00000002 & R1)!=0);  \
                    N_local[2]+=(uint##ISIZE##_t)(( 0x00000004 & R1)!=0);  \
                    N_local[3]+=(uint##ISIZE##_t)(( 0x00000008 & R1)!=0);  \
                    N_local[4]+=(uint##ISIZE##_t)(( 0x00000010 & R1)!=0);  \
                    N_local[5]+=(uint##ISIZE##_t)(( 0x00000020 & R1)!=0);  \
                    N_local[6]+=(uint##ISIZE##_t)(( 0x00000040 & R1)!=0);  \
                    N_local[7]+=(uint##ISIZE##_t)(( 0x00000080 & R1)!=0);  \
                    N_local[8]+=(uint##ISIZE##_t)(( 0x00000100 & R1)!=0);  \
                    N_local[9]+=(uint##ISIZE##_t)(( 0x00000200 & R1)!=0);  \
                    N_local[10]+=(uint##ISIZE##_t)(( 0x00000400 & R1)!=0);  \
                    N_local[11]+=(uint##ISIZE##_t)(( 0x00000800 & R1)!=0);  \
                    N_local[12]+=(uint##ISIZE##_t)(( 0x00001000 & R1)!=0);  \
                    N_local[13]+=(uint##ISIZE##_t)(( 0x00002000 & R1)!=0);  \
                    N_local[14]+=(uint##ISIZE##_t)(( 0x00004000 & R1)!=0);  \
                    N_local[15]+=(uint##ISIZE##_t)(( 0x00008000 & R1)!=0);  \
                    N_local[16]+=(uint##ISIZE##_t)(( 0x00010000 & R1)!=0);  \
                    N_local[17]+=(uint##ISIZE##_t)(( 0x00020000 & R1)!=0);  \
                    N_local[18]+=(uint##ISIZE##_t)(( 0x00040000 & R1)!=0);  \
                    N_local[19]+=(uint##ISIZE##_t)(( 0x00080000 & R1)!=0);  \
                    N_local[20]+=(uint##ISIZE##_t)(( 0x00100000 & R1)!=0);  \
                    N_local[21]+=(uint##ISIZE##_t)(( 0x00200000 & R1)!=0);  \
                    N_local[22]+=(uint##ISIZE##_t)(( 0x00400000 & R1)!=0);  \
                    N_local[23]+=(uint##ISIZE##_t)(( 0x00800000 & R1)!=0);  \
                    N_local[24]+=(uint##ISIZE##_t)(( 0x01000000 & R1)!=0);  \
                    N_local[25]+=(uint##ISIZE##_t)(( 0x02000000 & R1)!=0);  \
                    N_local[26]+=(uint##ISIZE##_t)(( 0x04000000 & R1)!=0);  \
                    N_local[27]+=(uint##ISIZE##_t)(( 0x08000000 & R1)!=0);  \
                    N_local[28]+=(uint##ISIZE##_t)(( 0x10000000 & R1)!=0);  \
                    N_local[29]+=(uint##ISIZE##_t)(( 0x20000000 & R1)!=0);  \
                    N_local[30]+=(uint##ISIZE##_t)(( 0x40000000 & R1)!=0);  \
                    N_local[31]+=(uint##ISIZE##_t)(( 0x80000000 & R1)!=0);  \
                    j2=32;                                                                                     \
                    for (i2=1;i2<m;i2++) {                                                                    \
                        R1=R[j1+i2];                                                                  \
                        N_local[j2+0]+=(uint##ISIZE##_t)(( 0x00000001 & R1)!=0);  \
                        N_local[j2+1]+=(uint##ISIZE##_t)(( 0x00000002 & R1)!=0);  \
                        N_local[j2+2]+=(uint##ISIZE##_t)(( 0x00000004 & R1)!=0);  \
                        N_local[j2+3]+=(uint##ISIZE##_t)(( 0x00000008 & R1)!=0);  \
                        N_local[j2+4]+=(uint##ISIZE##_t)(( 0x00000010 & R1)!=0);  \
                        N_local[j2+5]+=(uint##ISIZE##_t)(( 0x00000020 & R1)!=0);  \
                        N_local[j2+6]+=(uint##ISIZE##_t)(( 0x00000040 & R1)!=0);  \
                        N_local[j2+7]+=(uint##ISIZE##_t)(( 0x00000080 & R1)!=0);  \
                        N_local[j2+8]+=(uint##ISIZE##_t)(( 0x00000100 & R1)!=0);  \
                        N_local[j2+9]+=(uint##ISIZE##_t)(( 0x00000200 & R1)!=0);  \
                        N_local[j2+10]+=(uint##ISIZE##_t)(( 0x00000400 & R1)!=0);  \
                        N_local[j2+11]+=(uint##ISIZE##_t)(( 0x00000800 & R1)!=0);  \
                        N_local[j2+12]+=(uint##ISIZE##_t)(( 0x00001000 & R1)!=0);  \
                        N_local[j2+13]+=(uint##ISIZE##_t)(( 0x00002000 & R1)!=0);  \
                        N_local[j2+14]+=(uint##ISIZE##_t)(( 0x00004000 & R1)!=0);  \
                        N_local[j2+15]+=(uint##ISIZE##_t)(( 0x00008000 & R1)!=0);  \
                        N_local[j2+16]+=(uint##ISIZE##_t)(( 0x00010000 & R1)!=0);  \
                        N_local[j2+17]+=(uint##ISIZE##_t)(( 0x00020000 & R1)!=0);  \
                        N_local[j2+18]+=(uint##ISIZE##_t)(( 0x00040000 & R1)!=0);  \
                        N_local[j2+19]+=(uint##ISIZE##_t)(( 0x00080000 & R1)!=0);  \
                        N_local[j2+20]+=(uint##ISIZE##_t)(( 0x00100000 & R1)!=0);  \
                        N_local[j2+21]+=(uint##ISIZE##_t)(( 0x00200000 & R1)!=0);  \
                        N_local[j2+22]+=(uint##ISIZE##_t)(( 0x00400000 & R1)!=0);  \
                        N_local[j2+23]+=(uint##ISIZE##_t)(( 0x00800000 & R1)!=0);  \
                        N_local[j2+24]+=(uint##ISIZE##_t)(( 0x01000000 & R1)!=0);  \
                        N_local[j2+25]+=(uint##ISIZE##_t)(( 0x02000000 & R1)!=0);  \
                        N_local[j2+26]+=(uint##ISIZE##_t)(( 0x04000000 & R1)!=0);  \
                        N_local[j2+27]+=(uint##ISIZE##_t)(( 0x08000000 & R1)!=0);  \
                        N_local[j2+28]+=(uint##ISIZE##_t)(( 0x10000000 & R1)!=0);  \
                        N_local[j2+29]+=(uint##ISIZE##_t)(( 0x20000000 & R1)!=0);  \
                        N_local[j2+30]+=(uint##ISIZE##_t)(( 0x40000000 & R1)!=0);  \
                        N_local[j2+31]+=(uint##ISIZE##_t)(( 0x80000000 & R1)!=0);  \
                        j2+=32;                                                                            \
                    }                                                                                         \
                }                                                                                             \
            }

# define FOR_LOOP_64(ISIZE)                                                                             \
            for (i1=0;i1<ny;i1++) {                                                                           \
                j1=(uint64_t) i_perm[i1];                                                             \
                j1*=m;                                                                                        \
                i3=0;                                                                                         \
                while ((( R_subset[i3]| R[j1+i3])==0xffffffffffffffff)&&(i3<m)) {i3++;}                       \
                if (i3==m) {                                                                                  \
                    R1=R[j1];                                                                  \
                    N_local[0]+=(uint##ISIZE##_t)(( 0x0000000000000001 & R1)!=0);  \
                    N_local[1]+=(uint##ISIZE##_t)(( 0x0000000000000002 & R1)!=0);  \
                    N_local[2]+=(uint##ISIZE##_t)(( 0x0000000000000004 & R1)!=0);  \
                    N_local[3]+=(uint##ISIZE##_t)(( 0x0000000000000008 & R1)!=0);  \
                    N_local[4]+=(uint##ISIZE##_t)(( 0x0000000000000010 & R1)!=0);  \
                    N_local[5]+=(uint##ISIZE##_t)(( 0x0000000000000020 & R1)!=0);  \
                    N_local[6]+=(uint##ISIZE##_t)(( 0x0000000000000040 & R1)!=0);  \
                    N_local[7]+=(uint##ISIZE##_t)(( 0x0000000000000080 & R1)!=0);  \
                    N_local[8]+=(uint##ISIZE##_t)(( 0x0000000000000100 & R1)!=0);  \
                    N_local[9]+=(uint##ISIZE##_t)(( 0x0000000000000200 & R1)!=0);  \
                    N_local[10]+=(uint##ISIZE##_t)(( 0x0000000000000400 & R1)!=0);  \
                    N_local[11]+=(uint##ISIZE##_t)(( 0x0000000000000800 & R1)!=0);  \
                    N_local[12]+=(uint##ISIZE##_t)(( 0x0000000000001000 & R1)!=0);  \
                    N_local[13]+=(uint##ISIZE##_t)(( 0x0000000000002000 & R1)!=0);  \
                    N_local[14]+=(uint##ISIZE##_t)(( 0x0000000000004000 & R1)!=0);  \
                    N_local[15]+=(uint##ISIZE##_t)(( 0x0000000000008000 & R1)!=0);  \
                    N_local[16]+=(uint##ISIZE##_t)(( 0x0000000000010000 & R1)!=0);  \
                    N_local[17]+=(uint##ISIZE##_t)(( 0x0000000000020000 & R1)!=0);  \
                    N_local[18]+=(uint##ISIZE##_t)(( 0x0000000000040000 & R1)!=0);  \
                    N_local[19]+=(uint##ISIZE##_t)(( 0x0000000000080000 & R1)!=0);  \
                    N_local[20]+=(uint##ISIZE##_t)(( 0x0000000000100000 & R1)!=0);  \
                    N_local[21]+=(uint##ISIZE##_t)(( 0x0000000000200000 & R1)!=0);  \
                    N_local[22]+=(uint##ISIZE##_t)(( 0x0000000000400000 & R1)!=0);  \
                    N_local[23]+=(uint##ISIZE##_t)(( 0x0000000000800000 & R1)!=0);  \
                    N_local[24]+=(uint##ISIZE##_t)(( 0x0000000001000000 & R1)!=0);  \
                    N_local[25]+=(uint##ISIZE##_t)(( 0x0000000002000000 & R1)!=0);  \
                    N_local[26]+=(uint##ISIZE##_t)(( 0x0000000004000000 & R1)!=0);  \
                    N_local[27]+=(uint##ISIZE##_t)(( 0x0000000008000000 & R1)!=0);  \
                    N_local[28]+=(uint##ISIZE##_t)(( 0x0000000010000000 & R1)!=0);  \
                    N_local[29]+=(uint##ISIZE##_t)(( 0x0000000020000000 & R1)!=0);  \
                    N_local[30]+=(uint##ISIZE##_t)(( 0x0000000040000000 & R1)!=0);  \
                    N_local[31]+=(uint##ISIZE##_t)(( 0x0000000080000000 & R1)!=0);  \
                    N_local[32]+=(uint##ISIZE##_t)(( 0x0000000100000000 & R1)!=0);  \
                    N_local[33]+=(uint##ISIZE##_t)(( 0x0000000200000000 & R1)!=0);  \
                    N_local[34]+=(uint##ISIZE##_t)(( 0x0000000400000000 & R1)!=0);  \
                    N_local[35]+=(uint##ISIZE##_t)(( 0x0000000800000000 & R1)!=0);  \
                    N_local[36]+=(uint##ISIZE##_t)(( 0x0000001000000000 & R1)!=0);  \
                    N_local[37]+=(uint##ISIZE##_t)(( 0x0000002000000000 & R1)!=0);  \
                    N_local[38]+=(uint##ISIZE##_t)(( 0x0000004000000000 & R1)!=0);  \
                    N_local[39]+=(uint##ISIZE##_t)(( 0x0000008000000000 & R1)!=0);  \
                    N_local[40]+=(uint##ISIZE##_t)(( 0x0000010000000000 & R1)!=0);  \
                    N_local[41]+=(uint##ISIZE##_t)(( 0x0000020000000000 & R1)!=0);  \
                    N_local[42]+=(uint##ISIZE##_t)(( 0x0000040000000000 & R1)!=0);  \
                    N_local[43]+=(uint##ISIZE##_t)(( 0x0000080000000000 & R1)!=0);  \
                    N_local[44]+=(uint##ISIZE##_t)(( 0x0000100000000000 & R1)!=0);  \
                    N_local[45]+=(uint##ISIZE##_t)(( 0x0000200000000000 & R1)!=0);  \
                    N_local[46]+=(uint##ISIZE##_t)(( 0x0000400000000000 & R1)!=0);  \
                    N_local[47]+=(uint##ISIZE##_t)(( 0x0000800000000000 & R1)!=0);  \
                    N_local[48]+=(uint##ISIZE##_t)(( 0x0001000000000000 & R1)!=0);  \
                    N_local[49]+=(uint##ISIZE##_t)(( 0x0002000000000000 & R1)!=0);  \
                    N_local[50]+=(uint##ISIZE##_t)(( 0x0004000000000000 & R1)!=0);  \
                    N_local[51]+=(uint##ISIZE##_t)(( 0x0008000000000000 & R1)!=0);  \
                    N_local[52]+=(uint##ISIZE##_t)(( 0x0010000000000000 & R1)!=0);  \
                    N_local[53]+=(uint##ISIZE##_t)(( 0x0020000000000000 & R1)!=0);  \
                    N_local[54]+=(uint##ISIZE##_t)(( 0x0040000000000000 & R1)!=0);  \
                    N_local[55]+=(uint##ISIZE##_t)(( 0x0080000000000000 & R1)!=0);  \
                    N_local[56]+=(uint##ISIZE##_t)(( 0x0100000000000000 & R1)!=0);  \
                    N_local[57]+=(uint##ISIZE##_t)(( 0x0200000000000000 & R1)!=0);  \
                    N_local[58]+=(uint##ISIZE##_t)(( 0x0400000000000000 & R1)!=0);  \
                    N_local[59]+=(uint##ISIZE##_t)(( 0x0800000000000000 & R1)!=0);  \
                    N_local[60]+=(uint##ISIZE##_t)(( 0x1000000000000000 & R1)!=0);  \
                    N_local[61]+=(uint##ISIZE##_t)(( 0x2000000000000000 & R1)!=0);  \
                    N_local[62]+=(uint##ISIZE##_t)(( 0x4000000000000000 & R1)!=0);  \
                    N_local[63]+=(uint##ISIZE##_t)(( 0x8000000000000000 & R1)!=0);  \
                    j2=64;                                                                                     \
                    for (i2=1;i2<m;i2++) {                                                                    \
                        R1=R[j1+i2];                                                                  \
                        N_local[j2+0]+=(uint##ISIZE##_t)(( 0x0000000000000001 & R1)!=0);  \
                        N_local[j2+1]+=(uint##ISIZE##_t)(( 0x0000000000000002 & R1)!=0);  \
                        N_local[j2+2]+=(uint##ISIZE##_t)(( 0x0000000000000004 & R1)!=0);  \
                        N_local[j2+3]+=(uint##ISIZE##_t)(( 0x0000000000000008 & R1)!=0);  \
                        N_local[j2+4]+=(uint##ISIZE##_t)(( 0x0000000000000010 & R1)!=0);  \
                        N_local[j2+5]+=(uint##ISIZE##_t)(( 0x0000000000000020 & R1)!=0);  \
                        N_local[j2+6]+=(uint##ISIZE##_t)(( 0x0000000000000040 & R1)!=0);  \
                        N_local[j2+7]+=(uint##ISIZE##_t)(( 0x0000000000000080 & R1)!=0);  \
                        N_local[j2+8]+=(uint##ISIZE##_t)(( 0x0000000000000100 & R1)!=0);  \
                        N_local[j2+9]+=(uint##ISIZE##_t)(( 0x0000000000000200 & R1)!=0);  \
                        N_local[j2+10]+=(uint##ISIZE##_t)(( 0x0000000000000400 & R1)!=0);  \
                        N_local[j2+11]+=(uint##ISIZE##_t)(( 0x0000000000000800 & R1)!=0);  \
                        N_local[j2+12]+=(uint##ISIZE##_t)(( 0x0000000000001000 & R1)!=0);  \
                        N_local[j2+13]+=(uint##ISIZE##_t)(( 0x0000000000002000 & R1)!=0);  \
                        N_local[j2+14]+=(uint##ISIZE##_t)(( 0x0000000000004000 & R1)!=0);  \
                        N_local[j2+15]+=(uint##ISIZE##_t)(( 0x0000000000008000 & R1)!=0);  \
                        N_local[j2+16]+=(uint##ISIZE##_t)(( 0x0000000000010000 & R1)!=0);  \
                        N_local[j2+17]+=(uint##ISIZE##_t)(( 0x0000000000020000 & R1)!=0);  \
                        N_local[j2+18]+=(uint##ISIZE##_t)(( 0x0000000000040000 & R1)!=0);  \
                        N_local[j2+19]+=(uint##ISIZE##_t)(( 0x0000000000080000 & R1)!=0);  \
                        N_local[j2+20]+=(uint##ISIZE##_t)(( 0x0000000000100000 & R1)!=0);  \
                        N_local[j2+21]+=(uint##ISIZE##_t)(( 0x0000000000200000 & R1)!=0);  \
                        N_local[j2+22]+=(uint##ISIZE##_t)(( 0x0000000000400000 & R1)!=0);  \
                        N_local[j2+23]+=(uint##ISIZE##_t)(( 0x0000000000800000 & R1)!=0);  \
                        N_local[j2+24]+=(uint##ISIZE##_t)(( 0x0000000001000000 & R1)!=0);  \
                        N_local[j2+25]+=(uint##ISIZE##_t)(( 0x0000000002000000 & R1)!=0);  \
                        N_local[j2+26]+=(uint##ISIZE##_t)(( 0x0000000004000000 & R1)!=0);  \
                        N_local[j2+27]+=(uint##ISIZE##_t)(( 0x0000000008000000 & R1)!=0);  \
                        N_local[j2+28]+=(uint##ISIZE##_t)(( 0x0000000010000000 & R1)!=0);  \
                        N_local[j2+29]+=(uint##ISIZE##_t)(( 0x0000000020000000 & R1)!=0);  \
                        N_local[j2+30]+=(uint##ISIZE##_t)(( 0x0000000040000000 & R1)!=0);  \
                        N_local[j2+31]+=(uint##ISIZE##_t)(( 0x0000000080000000 & R1)!=0);  \
                        N_local[j2+32]+=(uint##ISIZE##_t)(( 0x0000000100000000 & R1)!=0);  \
                        N_local[j2+33]+=(uint##ISIZE##_t)(( 0x0000000200000000 & R1)!=0);  \
                        N_local[j2+34]+=(uint##ISIZE##_t)(( 0x0000000400000000 & R1)!=0);  \
                        N_local[j2+35]+=(uint##ISIZE##_t)(( 0x0000000800000000 & R1)!=0);  \
                        N_local[j2+36]+=(uint##ISIZE##_t)(( 0x0000001000000000 & R1)!=0);  \
                        N_local[j2+37]+=(uint##ISIZE##_t)(( 0x0000002000000000 & R1)!=0);  \
                        N_local[j2+38]+=(uint##ISIZE##_t)(( 0x0000004000000000 & R1)!=0);  \
                        N_local[j2+39]+=(uint##ISIZE##_t)(( 0x0000008000000000 & R1)!=0);  \
                        N_local[j2+40]+=(uint##ISIZE##_t)(( 0x0000010000000000 & R1)!=0);  \
                        N_local[j2+41]+=(uint##ISIZE##_t)(( 0x0000020000000000 & R1)!=0);  \
                        N_local[j2+42]+=(uint##ISIZE##_t)(( 0x0000040000000000 & R1)!=0);  \
                        N_local[j2+43]+=(uint##ISIZE##_t)(( 0x0000080000000000 & R1)!=0);  \
                        N_local[j2+44]+=(uint##ISIZE##_t)(( 0x0000100000000000 & R1)!=0);  \
                        N_local[j2+45]+=(uint##ISIZE##_t)(( 0x0000200000000000 & R1)!=0);  \
                        N_local[j2+46]+=(uint##ISIZE##_t)(( 0x0000400000000000 & R1)!=0);  \
                        N_local[j2+47]+=(uint##ISIZE##_t)(( 0x0000800000000000 & R1)!=0);  \
                        N_local[j2+48]+=(uint##ISIZE##_t)(( 0x0001000000000000 & R1)!=0);  \
                        N_local[j2+49]+=(uint##ISIZE##_t)(( 0x0002000000000000 & R1)!=0);  \
                        N_local[j2+50]+=(uint##ISIZE##_t)(( 0x0004000000000000 & R1)!=0);  \
                        N_local[j2+51]+=(uint##ISIZE##_t)(( 0x0008000000000000 & R1)!=0);  \
                        N_local[j2+52]+=(uint##ISIZE##_t)(( 0x0010000000000000 & R1)!=0);  \
                        N_local[j2+53]+=(uint##ISIZE##_t)(( 0x0020000000000000 & R1)!=0);  \
                        N_local[j2+54]+=(uint##ISIZE##_t)(( 0x0040000000000000 & R1)!=0);  \
                        N_local[j2+55]+=(uint##ISIZE##_t)(( 0x0080000000000000 & R1)!=0);  \
                        N_local[j2+56]+=(uint##ISIZE##_t)(( 0x0100000000000000 & R1)!=0);  \
                        N_local[j2+57]+=(uint##ISIZE##_t)(( 0x0200000000000000 & R1)!=0);  \
                        N_local[j2+58]+=(uint##ISIZE##_t)(( 0x0400000000000000 & R1)!=0);  \
                        N_local[j2+59]+=(uint##ISIZE##_t)(( 0x0800000000000000 & R1)!=0);  \
                        N_local[j2+60]+=(uint##ISIZE##_t)(( 0x1000000000000000 & R1)!=0);  \
                        N_local[j2+61]+=(uint##ISIZE##_t)(( 0x2000000000000000 & R1)!=0);  \
                        N_local[j2+62]+=(uint##ISIZE##_t)(( 0x4000000000000000 & R1)!=0);  \
                        N_local[j2+63]+=(uint##ISIZE##_t)(( 0x8000000000000000 & R1)!=0);  \
                        j2+=64;                                                                            \
                    }                                                                                         \
                }                                                                                             \
            }

# define COPY_LOCAL_VARS                        \
            {                                   \
                for (i2=0;i2<m*64;i2++) {       \
                    N[i2]+=(double)N_local[i2]; \
                }                               \
                free(N_local);                  \
            }

            
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int64_t ny;
    double *N;
    uint32_t i;
    bool    mode_64bit;

    
    if (nrhs != 4) {
        mexErrMsgTxt("Four input arguments required . . .");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("One output argument required . . .");
    }
    
    if (mxGetM(prhs[0])!=mxGetM(prhs[1])) {
        mexErrMsgTxt("Input arguments R_subset and R must have the same number of rows . . .");
    }

    ny= mxGetN(prhs[2]);
    ny*= mxGetM(prhs[2]);
    
    mode_64bit=(mxGetM(prhs[0])&1)==0;
    

    plhs[0]=mxCreateNumericMatrix(mxGetM(prhs[0])*64, 1, mxDOUBLE_CLASS, mxREAL);
    N = (double*) mxGetPr(plhs[0]);
    
    for (i=0;i<mxGetM(prhs[0])*64;i++) {N[i]=0;}
    
    
    uint16_t* max_threads = (uint16_t*) mxGetPr(prhs[3]);
    
    if (max_threads[0]<omp_get_num_procs()) {
        omp_set_num_threads(max_threads[0]);
    } else {
        omp_set_num_threads(omp_get_num_procs());        
    }
    
    mxClassID  cl = mxGetClassID(prhs[2]);
    if (cl==mxUINT32_CLASS) {
        if (mode_64bit) {
            VAR_INIT0(64,uint32_t)
            #pragma omp parallel  shared(N,ny) 
            {
                VAR_INIT(64,32)
                #pragma omp for
                FOR_LOOP_64(32)
                # pragma omp critical
                COPY_LOCAL_VARS
            }
        } else {
            VAR_INIT0(32,uint32_t)
            #pragma omp parallel  shared(N,ny) 
            {
                VAR_INIT(32,32)
                #pragma omp for
                FOR_LOOP_32(32)
                # pragma omp critical
                COPY_LOCAL_VARS
            }
        }
    } else {
        if (cl==mxDOUBLE_CLASS) {
            if (mode_64bit) {
                VAR_INIT0(64,double)
                #pragma omp parallel  shared(N,ny) 
                {
                    VAR_INIT(64,64)
                    #pragma omp for
                    FOR_LOOP_64(64)
                    # pragma omp critical
                    COPY_LOCAL_VARS
                }    
            } else {
                VAR_INIT0(32,double)
                #pragma omp parallel  shared(N,ny) 
                {
                    VAR_INIT(32,64)
                    #pragma omp for
                    FOR_LOOP_32(64)
                    # pragma omp critical
                    COPY_LOCAL_VARS
                }    
            }
        } else {
            mexErrMsgTxt("i_perm should be of type 'uint32' or 'double' . . .");
        }
    }
    
    
    return;
}