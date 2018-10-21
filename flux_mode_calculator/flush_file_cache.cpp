/*

This mex function is part of the MATLAB function CALCULATE_FLUX_MODES

flush_file_cache(filename)

'filename' will be flushed from the operating system's file cache; 
under windows also the associated memory will be freed

*/


#include "mex.h"

#ifdef _MSC_VER
#   include "stdio.h"
#   include "windows.h"
#elif __GNUC__
#   include <unistd.h>
#   include <sys/stat.h>
#   include <fcntl.h>
#endif






void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    char *x;

    
    if (nrhs != 1) {
        mexErrMsgTxt("One input argument expected . . .");
    }
    if (nlhs != 0) {
        mexErrMsgTxt("No output arguments allowed . . .");
    }
    if (mxIsChar(prhs[0]) != 1) {
        mexErrMsgTxt("Input must be a string . . .");
    }
    
 
    int n = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

    x=(char*) mxCalloc(n, sizeof(char));
	mxGetString(prhs[0],x,n);

#ifdef _MSC_VER  
    HANDLE hFile = CreateFile(x, GENERIC_WRITE, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile != INVALID_HANDLE_VALUE) {
        FlushFileBuffers(hFile);
        CloseHandle(hFile); 
    }
#elif __GNUC__
    int fd = open(x, O_WRONLY, 0660);
    if (fd != -1) {
        fsync (fd);
        close (fd);
    }
#endif
    
    return;
}




