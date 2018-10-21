/*

This mex function is part of the MATLAB function CALCULATE_FLUX_MODES

i = adjacency_test_pattern_tree_C(R,j_comb_cand,R_tree,i_node,i_leaf,max_threads);

with as input 
    R                       - the binary (non-zero) flux modes
    j_comb_cand             - indices of candidate adjacent flux mode combinations (after pruning)
    R_tree, i_node, i_leaf  - binary pattern tree of R
    max_threads             - maximum number of threads to use

and as output
    i                       - logical array with indices of truly adjacent flux modes


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

#define VAR_INIT0(RSIZE,JSIZE)                                                        \
        const uint32_t m=(uint32_t)mxGetM(prhs[0])*32/RSIZE;                    \
        uint##RSIZE##_t* const R = (uint##RSIZE##_t*) mxGetPr(prhs[0]);         \
        uint##RSIZE##_t* const R_tree = (uint##RSIZE##_t*) mxGetPr(prhs[2]);  \
        JSIZE* const j_comb = (JSIZE*) mxGetPr(prhs[1]);

        
#define VAR_INIT(RSIZE)                 \
        int64_t i_comb;                 \
        uint64_t j1,j2,i1,nm;           \
        uint32_t i,i2,i_treeX;          \
        uint32_t i_tree,n_stack;        \
        uint32_t i_tree_stack[1000];    \
        uint64_t j_tree;                \
        uint8_t a;                      \
        uint##RSIZE##_t *R_comb = (uint##RSIZE##_t*) malloc(m*sizeof(uint##RSIZE##_t));

//
// adjacency test is faster if ordered (right-branch) leafs are checked before unordered (left-branch) leafs
//
# define FOR_LOOP                                                                 \
        for (i_comb=0;i_comb<n_comb;i_comb++) {                                                 \
            j1=(uint64_t) j_comb[2*i_comb];                                                     \
            j1*=m;                                                                              \
            j2=(uint64_t) j_comb[2*i_comb+1];                                                   \
            j2*=m;                                                                              \
            for (i=0;i<m;i++) {                                                                 \
                R_comb[i]=~(R[j1+i]|R[j2+i]);                                   \
            }                                                                                   \
            a=0;                                                                                \
            i_tree=0;                                                                           \
            n_stack=0;                                                                          \
            i_tree_stack[0]=0;                                                                  \
            while ((a==0)&&(i_tree<n_tree)) {                                                   \
                if (i_tree<n_tree) {                                                            \
                    i=0;                                                                        \
                    j_tree=(uint64_t) i_tree;                                                   \
                    j_tree*=m;                                                                  \
                    while ((i<m)&&((R_comb[i] & R_tree[j_tree+i])==0)) {i++;};          \
                    if (i==m) {                                                                 \
                        n_stack++;                                                              \
                        i_tree_stack[n_stack]=i_tree;                                           \
                        i_tree++;                                                               \
                    } else {                                                                    \
                        i_tree=i_node[i_tree];                                                  \
                    }                                                                           \
                    i_treeX=i_tree_stack[n_stack];                                              \
                    while ((n_stack>0)&&(i_node[i_treeX]<=i_tree)&&(a==0)) {                    \
                        i1=i_leaf[i_treeX*2];                                                   \
                        i1*=m;                                                                  \
                        nm=i_leaf[i_treeX*2+1];                                                 \
                        nm*=m;                                                                  \
                        while ((a==0)&&(i1<=nm)) {                                              \
                            i2=0;                                                               \
                            while ((i2<m)&&((R_comb[i2] & R[i1+i2])==0)) {i2++;};       \
                            a=(i2==m)&&(i1!=j1)&&(i1!=j2);                                      \
                            i1+=m;                                                              \
                        }                                                                       \
                        n_stack--;                                                              \
                        i_treeX=i_tree_stack[n_stack];                                          \
                    }                                                                           \
                }                                                                               \
            }                                                                                   \
            i_comb_out[i_comb]=!a;                                                              \
        }                                                                                       \
        free(R_comb);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    uint64_t n_comb;
    uint32_t *i_node;
    int64_t *i_leaf;
    uint8_t *i_comb_out;
    
    bool    mode_64bit;

    
    if (nrhs != 6) {
        mexErrMsgTxt("Six input arguments required . . .");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("One output argument required . . .");
    }
    
    if (mxGetM(prhs[1])!=2) {
        mexErrMsgTxt("Candidate combinations must be a 2 by n matrix . . .");
    }
    if ((mxGetN(prhs[2])!=mxGetN(prhs[3]))||(mxGetN(prhs[2])!=mxGetN(prhs[4]))) {
        mexErrMsgTxt("R_tree, i_node and i_leaf must have the same number of rows . . .");
    }
    if (mxGetM(prhs[0])!=mxGetM(prhs[2])) {
        mexErrMsgTxt("R and R_tree must have the same number of rows . . .");
    }
    if (mxGetM(prhs[3])!=1) {
        mexErrMsgTxt("i_node must have one row . . .");
    }
    if (mxGetM(prhs[4])!=2) {
        mexErrMsgTxt("i_leaf must have two rows . . .");
    }

    mode_64bit=(mxGetM(prhs[0])&1)==0;
    
    n_comb=mxGetN(prhs[1]);
    
    i_node = (uint32_t*) mxGetPr(prhs[3]);
    i_leaf = (int64_t*) mxGetPr(prhs[4]);
    
    uint16_t* max_threads = (uint16_t*) mxGetPr(prhs[5]);
    
    if (max_threads[0]<omp_get_num_procs()) {
        omp_set_num_threads(max_threads[0]);
    } else {
        omp_set_num_threads(omp_get_num_procs());        
    }    
    
    const uint64_t n_tree = mxGetN(prhs[2]);
    

    plhs[0]=mxCreateNumericMatrix(n_comb, 1, mxUINT8_CLASS, mxREAL);
    i_comb_out = (uint8_t*) mxGetPr(plhs[0]);
    
    
    mxClassID  cl = mxGetClassID(prhs[1]);
    if (cl==mxUINT32_CLASS) {
        if (mode_64bit) {
            VAR_INIT0(64,uint32_t)
            #pragma omp parallel  shared(n_comb,i_comb_out,i_node,i_leaf) 
            {
                VAR_INIT(64)
                #pragma omp for schedule(dynamic,50)     
                FOR_LOOP
            }
        } else {
            VAR_INIT0(32,uint32_t)
            #pragma omp parallel  shared(n_comb,i_comb_out,i_node,i_leaf) 
            {
                VAR_INIT(32)
                #pragma omp for schedule(dynamic,50)     
                FOR_LOOP
            }            
        }
    } else {
        if (cl==mxDOUBLE_CLASS) {
            if (mode_64bit) {
                VAR_INIT0(64,double)
                #pragma omp parallel  shared(n_comb,i_comb_out,i_node,i_leaf) 
                {
                    VAR_INIT(64)
                    #pragma omp for schedule(dynamic,50)     
                    FOR_LOOP
                }
            } else {
                VAR_INIT0(32,double)
                #pragma omp parallel  shared(n_comb,i_comb_out,i_node,i_leaf) 
                {
                    VAR_INIT(32)
                    #pragma omp for schedule(dynamic,50)     
                    FOR_LOOP 
                }            
            }
        } else {
            mexErrMsgTxt("j_comb should be of type 'uint32' or 'double' . . .");
        }
    }
    return;
}