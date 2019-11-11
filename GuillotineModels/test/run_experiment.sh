#!/bin/bash

dtl="1.8E+6" # estimative of about 1 hour in the tiamats
instances=`ls ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/STS4{,s} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} P1_100_200_25_{1,2,3,4,5}`

for i in 1 2 3 4 5; do
	julia run_model.jl ---faithful2furini --no-{cut-position,redundant-cut} -seed "$i" --det-time-limit "$dtl" $instances | tee experiment_0_${i}_`date -Iminutes`.log
	julia run_model.jl --faithful2furini2016 --seed "$i" --det-time-limit "$dtl" $instances | tee experiment_1_${i}_`date -Iminutes`.log
	julia run_model.jl --faithful2furini2016 --round2disc --seed "$i" --det-time-limit "$dtl" $instances | tee experiment_2_${i}_`date -Iminutes`.log
	julia run_model.jl --seed "$i" --det-time-limit "$dtl" $instances | tee experiment_3__${i}_`date -Iminutes`.log
	julia run_model.jl --round2disc --seed "$i" --det-time-limit "$dtl" $instances | tee experiment_4_${i}_`date -Iminutes`.log
done

