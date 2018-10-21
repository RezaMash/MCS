#!/bin/bash
matlab -nodisplay -nodesktop -r "addpath('./flux_mode_calculator/'); model_name = 'e_coli_core'; sbml=true; MCS_by_null"
