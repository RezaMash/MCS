/*

This mex function is part of the MATLAB function ROUND_FRACTION

[A,B] = calculate_fraction(A_in,tol);

with as input 
    A_in                    - floating point matrix 
    tol                     - tolerance

and as output
    A,B                     - integer matrices such that A_in=A./B


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

#define max_number_continued_fraction 50





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double *A,*B,tol,x,a,b,h,k,h0,k0,e;
    long n,m,i,j;


    
    if (nrhs != 2) {
        mexErrMsgTxt("Two input arguments required . . .");
    }
    if (nlhs != 2) {
        mexErrMsgTxt("Two output arguments required . . .");
    }

    n=mxGetN(prhs[0])*mxGetM(prhs[0]);
    
        A=mxGetPr(prhs[1]);
        tol=A[0];
        plhs[1]=mxDuplicateArray(prhs[0]);
        if (mxIsEmpty(plhs[1])) {
            return;
        }
        B=mxGetPr(plhs[1]);
    
        plhs[0]=mxDuplicateArray(prhs[0]);
        if (mxIsEmpty(plhs[0])) {
            return;
        }
        A=mxGetPr(plhs[0]); 
        
        for (j=0;j<n;j++) {
            h=1;
            h0=0;
            k=0;
            k0=1;            
            x=A[j];
            
            a=floor(x);
            x-=a;
            b=h;
            h=h*a+h0;
            h0=b;
            b=k;
            k=k*a+k0;
            k0=b;
            
            if (x>tol) {
                i=1;
                e=tol+1;
                while ((e>tol)&&(i<max_number_continued_fraction)) {
                    x=1/x;
                    a=floor(x);
                    x-=a;
                    b=h;
                    h=h*a+h0;
                    h0=b;
                    b=k;
                    k=k*a+k0;
                    k0=b;
                    i++;
                    e=h/k-A[j];
                    if (e<0) {
                        e=-e;
                    }
                }
                if (x>=0.5) {
                    h+=h0;
                    k+=k0;
                }
                A[j]=h;
                B[j]=k;
            } else {
                A[j]=a;
                B[j]=1;
            }
        }



    
    return;
}