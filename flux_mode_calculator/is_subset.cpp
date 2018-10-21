/*

This mex function is part of the MATLAB function CALCULATE_FLUX_MODES

I = is_subset(R_subset,R,i_perm,max_threads);

with as input 
    R_subset                - bit patterns 
    R                       - the binary (non-zero) flux modes
    i_perm                  - indices of R over which to check the bit patterns
    max_threads             - maximum number of threads to use

and as output
    I                       - logical array with ones if R_subset is a subset of R(:,i_perm)


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

#define VAR_INIT            \
        uint32_t i2,i3,j2;  \
        int64_t i1;         \
        uint64_t j1a,j1;


# define FOR_LOOP                                                                 \
        for (i1=0;i1<ny;i1++) {                                                                 \
            j1=(uint64_t) i_perm[i1];                                                           \
            j1*=m;                                                                              \
            j1a=i1;                                                                             \
            j1a*=nx;                                                                            \
            for (i2=0;i2<nx;i2++) {                                                             \
                j2=i2*m;                                                                        \
                i3=0;                                                                           \
                while (((R_subset[i3+j2]&(~R[j1+i3]))==0)&&(i3<m)) {i3++;}      \
                I[j1a+i2]=i3==m;                                                                \
            }                                                                                   \
        }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    uint64_t nx,ny;
    bool    mode_64bit;
    uint8_t *I;


    
    if (nrhs != 4) {
        mexErrMsgTxt("Four input arguments required . . .");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("One output argument required . . .");
    }
    
    if (mxGetM(prhs[0])!=mxGetM(prhs[1])) {
        mexErrMsgTxt("Input arguments R_subset and R must have the same number of rows . . .");
    }

    nx=mxGetN(prhs[0]);
    ny=mxGetN(prhs[2])*mxGetM(prhs[2]);
    
    mode_64bit=(mxGetM(prhs[0])&1)==0;
        


    plhs[0]=mxCreateNumericMatrix(nx, ny, mxUINT8_CLASS, mxREAL);
    I = (uint8_t*) mxGetPr(plhs[0]);
    
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
            #pragma omp parallel  shared(I,nx,ny) 
            {
                VAR_INIT
                #pragma omp for
                FOR_LOOP
            }
        } else {
            VAR_INIT0(32,uint32_t)
            #pragma omp parallel  shared(I,nx,ny) 
            {
                VAR_INIT
                #pragma omp for
                FOR_LOOP
            }
        }
    } else {
        if (cl==mxDOUBLE_CLASS) {
            if (mode_64bit) {
                VAR_INIT0(64,double)
                #pragma omp parallel  shared(I,nx,ny) 
                {
                    VAR_INIT
                    #pragma omp for
                    FOR_LOOP
                }    
            } else {
                VAR_INIT0(32,double)
                #pragma omp parallel  shared(I,nx,ny) 
                {
                    VAR_INIT
                    #pragma omp for
                    FOR_LOOP
                }    
            }
        } else {
            mexErrMsgTxt("i_perm should be of type 'uint32' or 'double' . . .");
        }
    }

    
    
    return;
}