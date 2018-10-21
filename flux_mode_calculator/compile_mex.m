

% In order to install the glpk LP solver for performing network compression:
%  1) get the glpkcc.cpp file from OptiToolbox and copy it to the current directory
%  2) install glpk-448 and run Build_GLPK_with_VC9.bat
%  3) give glpk_dir the full path name of the directory in which glpk-4.48 is stored
%  4) eval(['mex -v -largeArrayDims glpkcc.cpp -I',glpk_dir,'\src -L',glpk_dir,'\w64 LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:libcmt.lib" -lglpk -lut -output glpk'])
%  5) rehash toolboxcache
%  6) run calculate_flux_modes with the option 'LPsolver' 'glpk'


% The mex functions called in the adjacency test make use of the OpenMP API
% to enable multithreading and should be compiled with the "/openmp" flag
% (MSVC) or the "-fopenmp" flag (GCC)

if ispc
    % windows platform, MSVC compiler
    mex COMPFLAGS="$COMPFLAGS /openmp"  -largeArrayDims -v count_subsets.cpp
    mex COMPFLAGS="$COMPFLAGS /openmp" -largeArrayDims -v is_subset.cpp
    mex COMPFLAGS="$COMPFLAGS /openmp"  -largeArrayDims -v adjacency_test_pattern_tree_C.cpp
    mex COMPFLAGS="$COMPFLAGS /openmp"  -largeArrayDims -v prune_candidates_rank_pattern_tree_C.cpp
    
    mex calculate_fraction.c        
    mex flush_file_cache.cpp
        
elseif isunix
    % unix platform, GCC compiler
    mex  CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims -v count_subsets.cpp
    mex CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims -v is_subset.cpp
    mex CXXFLAGS="\$CXXFLAGS -mpopcnt -fopenmp" LDFLAGS="\$LDFLAGS -mpopcnt -fopenmp" -largeArrayDims -v prune_candidates_rank_pattern_tree_C.cpp
    mex  CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims -v adjacency_test_pattern_tree_C.cpp
    
    mex calculate_fraction.c
    mex flush_file_cache.cpp
    
else
    error('Unknown operating system . . .')
end;

        