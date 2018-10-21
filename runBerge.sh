#!/bin/bash
if [ ! -d models/$1 ] ; then
  echo "Model $1 doesn't exist" >&2; exit 1
fi
if [ $3 = 0 ] ; then
  reactions=$(head -n 1 models/$1/full_matrix_human_readable.txt)
else
  reactions=$(head -n 1 models/$1/reduced_matrix_human_readable.txt)
fi

if [ $2 = -1 ] ; then
  log_name=$1_all-reactions_$3_berge-method_`date '+%Y-%m-%d_%H:%M:%S'`.csv
  printf "Procedure, Time in seconds, Number of MCSs\n" >> "logs/$log_name"
  matlab -nodisplay -nodesktop -r "reduced = $3; target_reaction = -1; model_name = '$1'; log_file = '$log_name'; sbml=true; BergeMCSsEnumator"
  printf "Total,=SUM(B2:B$(($reactions+2)))\nAverage,=AVERAGE(B2:B$(($reactions+2)))\n" >> "logs/$log_name"
else
  log_name=$1_$2_$3_berge-method_`date '+%Y-%m-%d_%H:%M:%S'`.csv
  printf "Procedure, Time in seconds, Number of MCSs\n" >> "logs/$log_name"
  matlab -nodisplay -nodesktop -r "reduced = $3; target_reaction = $2; model_name = '$1'; log_file = '$log_name'; sbml=true; BergeMCSsEnumator"
fi
