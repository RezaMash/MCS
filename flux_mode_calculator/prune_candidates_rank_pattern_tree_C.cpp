/*

This mex function is part of the MATLAB function CALCULATE_FLUX_MODES

[j_comb_cand,j1_inc,n_comb_cand] = prune_candidates_rank_pattern_tree_C(R,R_tree,i_node,i_leaf3,i_leaf1,th,Hamming_weights,max_blocks,j1_inc_in,filename,SSE4,max_threads);

with as input 
    R                       - the binary (non-zero) flux modes
    R_tree, i_node, i_leafX - binary pattern tree of R, with leafs for positive and negative flux modes
    th                      - pruning threshold
    Hamming_weights         - precomputed Hamming weights for numbers 0 .. 65535
    max_blocks              - maximum number of memory blocks of size block_size*4 that is allocated
    j1_inc_in               - index of positive flux modes that have been considered for pruning and adjacency test 
    filename                - location to which to write j_comb_cand if not enough free memory is available
    SSE4                    - 1 if internal popcnt instruction should be used - if supported by the processor; 0 otherwise
    max_threads             - maximum number of threads to use

and as output
    j_comb_cand             - indices of candidate adjacent flux mode combinations
    j1_inc                  - index of positive flux modes that have been considered for pruning and adjacency test
    n_comb_cand             - number of candidate adjacent flux mode combinations 


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



#ifdef _MSC_VER
#  include "intrin.h"
#elif __GNUC__
#  include "cpuid.h"
#endif


#define VAR_INIT0(RSIZE,JTYPE,JTYPE2)                                           \
        JTYPE2 *X[524288];                                                       \
        const uint32_t m=(uint32_t)mxGetM(prhs[0])*32/RSIZE;                    \
        uint##RSIZE##_t* const R = (uint##RSIZE##_t*) mxGetPr(prhs[0]);         \
        uint##RSIZE##_t* const R_tree = (uint##RSIZE##_t*) mxGetPr(prhs[1]);    \
        JTYPE *i_perm = (JTYPE*) malloc(n_tree*sizeof(JTYPE));                  \
        ii=0;                                                                   \
        nz_filt=0;                                                              \
        n_leaf_max=0;                                                           \
        for (i=0;i<n_tree;i++) {                                                \
            filt_tree[i]=1;                                                     \
            for (j=0;j<((i_leaf1[i*2+1]+1)-i_leaf1[i*2]);j++) {                 \
                filt_tree[i]=filt_tree[i]&(filt_in[ii+j]>=n_tree);              \
            }                                                                   \
            if (((i_leaf1[i*2+1]+1)-i_leaf1[i*2])>n_leaf_max) {                 \
                n_leaf_max=(uint32_t) ((i_leaf1[i*2+1]+1)-i_leaf1[i*2]);                   \
            }                                                                   \
            i_perm[i]=(JTYPE) ii;                                                       \
            if (filt_tree[i]==0) {nz_filt++;};                                  \
            ii+=(i_leaf1[i*2+1]+1)-i_leaf1[i*2];                                \
        }                                                                       \
        if (n1 != ii) {                                                         \
            free(i_perm);                                                       \
            free(filt_tree);                                                    \
            mexErrMsgTxt("j1_inc does not have the right length . . .");        \
        }
        

# define VAR_INIT(JTYPE)                                                               \
            uint32_t k=0,nn;                                                    \
            uint32_t n_write,i,n_leaf;                                          \
            int64_t j1,j3,j_tree1,j_tree3,i1,i3;                                   \
            int32_t  i_tree1,i_tree3_min,i_filt,x,i_leaf;                 \
            uint16_t nbits;                                                     \
            clock_t  t1,t2,dt;                                                  \
            dt=CLOCKS_PER_SEC*0.5;                                              \
            JTYPE *X_local = (JTYPE*) malloc(block_size*sizeof(JTYPE));\
            uint32_t *i_tree3 = (uint32_t*) malloc(n_leaf_max*sizeof(uint32_t));\
            i_tree1=-1;                                                         \
            x=0;
            
# define FOR_LOOP1(POPCNT,JTYPE) \
            for (i_filt=1;i_filt<=nz_filt;i_filt++) {                           \
                if (x<i_filt) {                                                 \
                    while (x<i_filt) {                                          \
                        i_tree1++;                                              \
                        x+=(filt_tree[i_tree1]==0);                             \
                    }                                                           \
                } else {                                                        \
                    while (x>i_filt) {                                          \
                        i_tree1--;                                              \
                        x-=(filt_tree[i_tree1]==0);                             \
                    }                                                           \
                }                                                               \
                j_tree1=i_tree1;                                                \
                j_tree1*=m;                                                     \
                n_leaf=(uint32_t) ((i_leaf1[i_tree1*2+1]+1)-i_leaf1[i_tree1*2]); \
                for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {                         \
                    i_tree3[i_leaf]=filt[i_perm[i_tree1]+i_leaf];               \
                }                                                               \
                i_tree3_min=n_tree;                                             \
                for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {                         \
                    if (i_tree3[i_leaf]<i_tree3_min) {                          \
                        i_tree3_min=i_tree3[i_leaf];                            \
                    }                                                           \
                }                                                               \
                while ((i_tree3_min<n_tree)&&(nX<=max_nblock)) {                \
                    j_tree3=i_tree3_min;                                        \
                    j_tree3*=m;                                                 \
                    nbits=0;                                                    \
                    for (i=0;i<m;i++) {nbits+=(uint16_t) POPCNT(R_tree[j_tree3+i]|R_tree[j_tree1+i]);} \
                    if (th>=nbits) {                                            \
                        for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {                 \
                            if (i_tree3[i_leaf]==i_tree3_min) {                 \
                                i1=i_leaf+i_leaf1[i_tree1*2];                   \
                                j1=i1;                                          \
                                j1*=m;                                          \
                                nbits=0;                                        \
                                for (i=0;i<m;i++) {nbits+=(uint16_t) POPCNT(R_tree[j_tree3+i]|R[j1+i]);} \
                                if (th<nbits) {                                 \
                                    i_tree3[i_leaf]=i_node[i_tree3[i_leaf]];    \
                                }                                               \
                            }                                                   \
                        }                                                       \
                        for (i3=i_leaf3[i_tree3_min*2];i3<=i_leaf3[i_tree3_min*2+1];i3++) { \
                            j3=i3;                                              \
                            j3*=m;                                              \
                            nbits=0;                                            \
                            for (i=0;i<m;i++) {nbits+=(uint16_t) POPCNT(R_tree[j_tree1+i]|R[j3+i]); } \
                            if (th>=nbits) {                                    \
                                for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {         \
                                    if (i_tree3[i_leaf]==i_tree3_min) {         \
                                        i1=i_leaf+i_leaf1[i_tree1*2];           \
                                        j1=i1;                                  \
                                        j1*=m;                                  \
                                        nbits=0;                                \
                                        for (i=0;i<m;i++) {nbits+=(uint16_t) POPCNT(R[j3+i]|R[j1+i]);} \
                                        if (th>=nbits) {                        \
                                            X_local[k]=(JTYPE) i1;    \
                                            X_local[k+1]=(JTYPE) i3;  \
                                            k+=2;                               \
                                            if (k<block_size) {                 \
                                            } else {


# define FOR_LOOP2(JTYPE) \
                                             {                                  \
                                                    X[nX] = (JTYPE*) malloc(k*sizeof(JTYPE));                         \
                                                    if (X[nX]==NULL) {                                                      \
                                                        pFile=fopen(j_comb_cand_filename,"ab");                             \
                                                        for (i=0;i<nX;i++) {                                                \
                                                            n_write=(uint32_t)fwrite(X[i],sizeof(JTYPE),mX[i],pFile);    \
                                                            nn=0;                                                           \
                                                            while ((n_write==0)&(nn<60)) {                                  \
                                                               t1=clock();                                                  \
                                                               t2=t1;                                                       \
                                                               while ((t2-t1)<dt) t2=clock();                               \
                                                               n_write=(uint32_t)fwrite(X[i],sizeof(JTYPE),mX[i],pFile); \
                                                               nn++;                                                        \
                                                            }                                                               \
                                                            free(X[i]);                                                     \
                                                        }                                                                   \
                                                        fclose(pFile);                                                      \
                                                        nX=0;                                                               \
                                                        write_to_disk=true;                                                 \
                                                        X[nX] = (JTYPE*) malloc(k*sizeof(JTYPE));                     \
                                                    }                           \
                                                    for (i=0;i<k;i++) {         \
                                                        X[nX][i]=(JTYPE) X_local[i];    \
                                                    }                           \
                                                    mX[nX]=k;                   \
                                                    n_comb+=k;                  \
                                                    nX++;                       \
                                                    k=0;                        \
                                                }                               \
                                            }                                   \
                                        }                                       \
                                    }                                           \
                                }                                               \
                            }                                                   \
                        }                                                       \
                        for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {                 \
                            if (i_tree3[i_leaf]==i_tree3_min) {                 \
                                i_tree3[i_leaf]++;                              \
                            }                                                   \
                        }                                                       \
                    } else {                                                    \
                        for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {                 \
                            if (i_tree3[i_leaf]==i_tree3_min) {                 \
                                i_tree3[i_leaf]=i_node[i_tree3[i_leaf]];        \
                            }                                                   \
                        }                                                       \
                    }                                                           \
                    i_tree3_min=n_tree;                                         \
                    for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {                     \
                        if (i_tree3[i_leaf]<i_tree3_min) {                      \
                            i_tree3_min=i_tree3[i_leaf];                        \
                        }                                                       \
                    }                                                           \
                }                                                               \
                for (i_leaf=0;i_leaf<n_leaf;i_leaf++) {                         \
                    filt[i_leaf+i_perm[i_tree1]]=i_tree3[i_leaf];               \
                }                                                               \
            }

# define COPY_LOCAL_VARS(JTYPE)                                                 \
            {                                                                   \
                X[nX] = (JTYPE*) malloc(k*sizeof(JTYPE));                 \
                if (X[nX]==NULL) {                                              \
                    pFile=fopen(j_comb_cand_filename,"ab");                     \
                    for (i=0;i<nX;i++) {                                        \
                        n_write=(uint32_t)fwrite(X[i],sizeof(JTYPE),mX[i],pFile); \
                        nn=0;                                                   \
                        while ((n_write==0)&(nn<60)) {                          \
                            t1=clock();                                         \
                            t2=t1;                                              \
                            while ((t2-t1)<dt) t2=clock();                      \
                            n_write=(uint32_t)fwrite(X[i],sizeof(JTYPE),mX[i],pFile);    \
                            nn++;                                               \
                        }                                                       \
                        free(X[i]);                                             \
                    }                                                           \
                    fclose(pFile);                                              \
                    nX=0;                                                       \
                    write_to_disk=true;                                         \
                    X[nX] = (JTYPE*) malloc(k*sizeof(JTYPE));             \
                }                                                               \
                for (i=0;i<k;i++) {                                             \
                    X[nX][i]=(JTYPE) X_local[i];                                        \
                }                                                               \
                mX[nX]=k;                                                       \
                n_comb+=k;                                                      \
                nX++;                                                           \
                free(X_local);                                                  \
                free(i_tree3);                                                  \
            }

#define POST_PROCESSING(JTYPE,JTYPE2)                                           \
    n_comb=n_comb/2;                                                            \
    if (!write_to_disk) {                                                       \
        JTYPE *dummy = (JTYPE*) malloc(2*n_comb*(sizeof(JTYPE)+1));    \
        if (dummy==NULL) {                                                      \
            write_to_disk=true;                                                 \
        } else {                                                                \
            free(dummy);                                                        \
        }                                                                       \
    }                                                                           \
    if (!write_to_disk) {                                                       \
        plhs[0]=mxCreateNumericMatrix(2, n_comb, mx##JTYPE2##_CLASS, mxREAL);       \
        JTYPE* j_comb = (JTYPE*) mxGetPr(plhs[0]);                        \
        jj=0;                                                                   \
        for (i=0;i<nX;i++) {                                                    \
            jj+=mX[i];                                                          \
        }                                                                       \
        if (jj!=2*n_comb) {                                                     \
            mexPrintf("\n");                                                    \
            mexPrintf("Number of combined rays: %d\n",n_comb*2);                \
            mexPrintf("Memory block sizes: \n");                                \
            jj=0;                                                               \
            for (i=0;i<nX;i++) {                                                \
                jj+=mX[i];                                                      \
                mexPrintf("%4d %11d %11d %11d\n",i+1,mX[i],jj,int64_t(n_comb*2)-int64_t(jj)); \
            }                                                                   \
            mexPrintf("\n");                                                    \
            mexErrMsgTxt("Total size of memory blocks does not correspond to n_comb . . ."); \
        }                                                                       \
        jj=0;                                                                   \
        for (i=0;i<nX;i++) {                                                    \
            for (j=0;j<mX[i];j++) {                                             \
                j_comb[jj+j]=X[i][j];                                           \
            }                                                                   \
            jj+=mX[i];                                                          \
            free(X[i]);                                                         \
        }                                                                       \
    } else {                                                                    \
        pFile=fopen(j_comb_cand_filename,"ab");                                 \
        for (j=0;j<nX;j++) {                                                    \
            fwrite(X[j],sizeof(JTYPE),mX[j],pFile);                          \
            free(X[j]);                                                         \
        }                                                                       \
        fclose(pFile);                                                          \
        plhs[0]=mxCreateNumericMatrix(1, 1, mx##JTYPE2##_CLASS, mxREAL);            \
        JTYPE* j_comb = (JTYPE*) mxGetPr(plhs[0]);                        \
    }                                                                           \
    free(i_perm);


#include "mex.h"
#include "stdio.h"
#include "math.h"
#include "omp.h"
#include "stdint.h"
#include "time.h"


#define block_size 2097156

uint16_t popcount_global[65536];

uint16_t popcount_lookup(uint16_t x) {
    return popcount_global[x];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    FILE *pFile;
    uint32_t nz_filt,i,j,nX,n_leaf_max;
    int64_t n_comb,jj,ii;
    uint32_t mX[524288];
    char *j_comb_cand_filename;
    bool write_to_disk;
    int CPUInfo[4] = {-1};
    bool mode_64bit;

    
    if (nrhs != 13) {
        mexErrMsgTxt("Thirteen input arguments required . . .");
    }
    if (nlhs != 3) {
        mexErrMsgTxt("Three output argument required . . .");
    }
    
    
    if (mxGetM(prhs[0])!=mxGetM(prhs[1])) {
        mexErrMsgTxt("R and R_tree must have the same number of rows . . .");
    }
    if ((mxGetN(prhs[1])!=mxGetN(prhs[2]))||(mxGetN(prhs[1])!=mxGetN(prhs[3]))||(mxGetN(prhs[1])!=mxGetN(prhs[4]))) {
        mexErrMsgTxt("R_tree, i_node and i_leaf must have the same number of columns . . .");
    }
    if (mxGetM(prhs[2])!=1) {
        mexErrMsgTxt("i_node must have one row . . .");
    }
    if ((mxGetM(prhs[3])!=2)||(mxGetM(prhs[4])!=2)) {
        mexErrMsgTxt("i_leaf must have two rows . . .");
    }
    if (mxGetM(prhs[6])*mxGetN(prhs[6])!=65536) {
        mexErrMsgTxt("The seventh argument must be a vector with the Hamming weights of all 16 bit integers . . .");
    }
    if (mxIsChar(prhs[9]) != 1) {
        mexErrMsgTxt("j_comb_cand_filename must be a string . . .");
    }

    uint32_t* p64_p = (uint32_t*) mxGetPr(prhs[12]);
    uint32_t p64 = p64_p[0];
    
    const int n_filename = (uint32_t) (mxGetM(prhs[9]) * mxGetN(prhs[9])) + 1;
    j_comb_cand_filename=(char*) mxCalloc(n_filename, sizeof(char));
	mxGetString(prhs[9],j_comb_cand_filename,n_filename);
    
    const uint32_t n_tree = (uint32_t) mxGetN(prhs[1]);
    mode_64bit=(mxGetM(prhs[0])&1)==0;
    
    uint32_t* const i_node = (uint32_t*) mxGetPr(prhs[2]);
    int64_t* const i_leaf3 = (int64_t*) mxGetPr(prhs[3]);
    int64_t* const i_leaf1 = (int64_t*) mxGetPr(prhs[4]);
  
    const uint64_t n1 = (mxGetM(prhs[8])*mxGetN(prhs[8]));
    
    plhs[1]=mxCreateNumericMatrix(1, (uint64_t)  n1, mxUINT32_CLASS, mxREAL);
    uint32_t* filt = (uint32_t*) mxGetPr(plhs[1]);
    uint32_t* filt_in = (uint32_t*) mxGetPr(prhs[8]);
    uint32_t *filt_tree = (uint32_t*) malloc(n_tree*sizeof(uint32_t));      
    for (ii=0;ii<n1;ii++) {                                                 
        filt[ii]=filt_in[ii];                                               
    }                                                                       

    uint32_t* SSE4_p = (uint32_t*) mxGetPr(prhs[10]);
    uint32_t SSE4 = SSE4_p[0];
    uint16_t* const th_p = (uint16_t*) mxGetPr(prhs[5]);
    const uint16_t th = th_p[0];
    uint8_t* const popcount = (uint8_t*) mxGetPr(prhs[6]);
    uint32_t* const max_nblock_p = (uint32_t*) mxGetPr(prhs[7]);
    
    uint16_t* max_threads_p = (uint16_t*) mxGetPr(prhs[11]);
    uint16_t max_threads = max_threads_p[0];
    
    if (max_threads>omp_get_num_procs()) {
        max_threads=omp_get_num_procs();        
    }
    
    omp_set_num_threads(max_threads);
    const uint32_t max_nblock = (max_nblock_p[0]-max_threads)/(p64+1);


    
    for (i=0;i<65536;i++) {popcount_global[i]=popcount[i];}

    for (i=0;i<65536;i++) {mX[i]=0;}
    nX=0;

    n_comb=0;


    write_to_disk=false;
    
    // check if popcnt is supported by the current processor
    if (SSE4!=0) {
        #ifdef _MSC_VER
            __cpuid(CPUInfo, 1);
        #elif __GNUC__
            __cpuid(1, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
        #endif
        SSE4=(CPUInfo[2] & (1 << 23));
    }
        
    if (SSE4==0) {
        if (p64==1) {
            VAR_INIT0(16,uint64_t,double)
            #pragma omp parallel  shared(n_comb,X,nX,mX,i_perm,n_leaf_max,filt_in,filt,write_to_disk,pFile,j_comb_cand_filename) 
            {
                VAR_INIT(uint64_t)
                #pragma omp for schedule(dynamic,1)
                FOR_LOOP1(popcount_lookup,uint64_t)
                # pragma omp critical
                FOR_LOOP2(double)
                # pragma omp critical
                COPY_LOCAL_VARS(double)
            }
            POST_PROCESSING(double,DOUBLE)
        } else {
            VAR_INIT0(16,uint32_t,uint32_t)
            #pragma omp parallel  shared(n_comb,X,nX,mX,i_perm,n_leaf_max,filt_in,filt,write_to_disk,pFile,j_comb_cand_filename) 
            {
                VAR_INIT(uint32_t)
                #pragma omp for schedule(dynamic,1)
                FOR_LOOP1(popcount_lookup,uint32_t)
                # pragma omp critical
                FOR_LOOP2(uint32_t)
                # pragma omp critical
                COPY_LOCAL_VARS(uint32_t)
            }            
            POST_PROCESSING(uint32_t,UINT32)
        }
    } else {
        
        if (mode_64bit) {
            if (p64==1) {
                VAR_INIT0(64,uint64_t,double)
                #pragma omp parallel  shared(n_comb,X,nX,mX,i_perm,n_leaf_max,filt_in,filt,write_to_disk,pFile,j_comb_cand_filename) 
                {
                    VAR_INIT(uint64_t)
                    #pragma omp for schedule(dynamic,1)
                    #ifdef _MSC_VER
                    FOR_LOOP1(__popcnt64,uint64_t)
                    #elif __GNUC__
                    FOR_LOOP1(__builtin_popcountll,uint64_t)
                    #endif
                    # pragma omp critical
                    FOR_LOOP2(double)
                    # pragma omp critical
                    COPY_LOCAL_VARS(double)
                }            
                POST_PROCESSING(double,DOUBLE)
            } else {
                VAR_INIT0(64,uint32_t,uint32_t)
                #pragma omp parallel  shared(n_comb,X,nX,mX,i_perm,n_leaf_max,filt_in,filt,write_to_disk,pFile,j_comb_cand_filename) 
                {
                    VAR_INIT(uint32_t)
                    #pragma omp for schedule(dynamic,1)
                    #ifdef _MSC_VER
                    FOR_LOOP1(__popcnt64,uint32_t)
                    #elif __GNUC__
                    FOR_LOOP1(__builtin_popcountll,uint32_t)
                    #endif
                    # pragma omp critical
                    FOR_LOOP2(uint32_t)
                    # pragma omp critical
                    COPY_LOCAL_VARS(uint32_t)
                }            
                POST_PROCESSING(uint32_t,UINT32)
            }
        } else {
            if (p64) {
                VAR_INIT0(32,uint64_t,double)
                #pragma omp parallel  shared(n_comb,X,nX,mX,i_perm,n_leaf_max,filt_in,filt,write_to_disk,pFile,j_comb_cand_filename) 
                {
                    VAR_INIT(uint64_t)
                    #pragma omp for schedule(dynamic,1)
                    #ifdef _MSC_VER
                    FOR_LOOP1(__popcnt,uint64_t)
                    #elif __GNUC__
                    FOR_LOOP1(__builtin_popcountl,uint64_t)
                    #endif
                    # pragma omp critical
                    FOR_LOOP2(double)
                    # pragma omp critical
                    COPY_LOCAL_VARS(double)
                }
                POST_PROCESSING(double,DOUBLE)
            } else {
                VAR_INIT0(32,uint32_t,uint32_t)
                #pragma omp parallel  shared(n_comb,X,nX,mX,i_perm,n_leaf_max,filt_in,filt,write_to_disk,pFile,j_comb_cand_filename) 
                {
                    VAR_INIT(uint32_t)
                    #pragma omp for schedule(dynamic,1)
                    #ifdef _MSC_VER
                    FOR_LOOP1(__popcnt,uint32_t)
                    #elif __GNUC__
                    FOR_LOOP1(__builtin_popcountl,uint32_t)
                    #endif
                    # pragma omp critical
                    FOR_LOOP2(uint32_t)
                    # pragma omp critical
                    COPY_LOCAL_VARS(uint32_t)
                }
                POST_PROCESSING(uint32_t,UINT32)
            }
        }
    }
    
    free(filt_tree);
  
    plhs[2]=mxCreateNumericMatrix(1,1, mxINT64_CLASS, mxREAL);    
    int64_t* n_comb1 = (int64_t*) mxGetPr(plhs[2]);
    n_comb1[0]=n_comb;
    
    return;
}