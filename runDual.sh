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
  log_name=$1_all-reactions_$3_dual-method_`date '+%Y-%m-%d_%H:%M:%S'`.csv
  printf "Target,EFMs finding time,post-processing time,Number of MCSs,total time in seconds\n" >> "logs/$log_name"
  for (( target=1; target<=$reactions; target++)); do
    matlab -nodisplay -nodesktop -r "addpath('./flux_mode_calculator/'); reduced = $3; target_reaction = $target; model_name = '$1'; log_file = '$log_name'; MCS_by_dual_Enhanced;"
    java -jar SuperSetRemoval.jar >> "logs/$log_name"
    printf ",=B$(($target+1))+C$(($target+1))\n" >> "logs/$log_name"
  done
  printf "Total,=SUM(B2:B$(($reactions+2))),=SUM(C2:C$(($reactions+2))),-,=SUM(E2:E$(($reactions+2)))\n" >> "logs/$log_name"
  printf "Average,=AVERAGE(B2:B$(($reactions+1))),=AVERAGE(C2:C$(($reactions+1))),-,=AVERAGE(E2:E$(($reactions+1)))\n" >> "logs/$log_name"
else
  log_name=$1_$2_$3_dual-method_`date '+%Y-%m-%d_%H:%M:%S'`.csv
  printf "Target,EFMs finding time,post-processing time,Number of MCSs,total time in seconds\n" >> "logs/$log_name"
  matlab -nodisplay -nodesktop -r "addpath('./flux_mode_calculator/'); reduced = $3; target_reaction = $2; model_name = '$1'; log_file = '$log_name'; MCS_by_dual_Enhanced;"
  java -jar SuperSetRemoval.jar >> "logs/$log_name"
  printf ",=B2+C2\n" >> "logs/$log_name"
fi
