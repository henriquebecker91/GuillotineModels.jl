#!/bin/bash

#dtl="1.8E+6" # estimative of about 1 hour in the tiamats
#dtl="3E+4" # estimative of about 1 minute in the tiamats
#tl="60.0"
tl="3600.0"
# ls is not used below because it reorders the instances by alphabetical order
instances=`echo ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/STS4{,s} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} ./instances/P1_100_200_25_{1,2,3,4,5}`
easy_for_furini=`echo ./instances/okp{1,4,5} ./instances/CU1 ./instances/STS4{,s} ./instances/gcut9`

julia run_model.jl --do-not-solve --faithful2furini2016 --no-{cut-position,redundant-cut} $instances 2>&1 | tee experiment_0_`date -Iminutes`.log
julia run_model.jl --do-not-solve --faithful2furini2016 --no-redundant-cut $instances 2>&1 | tee experiment_1_`date -Iminutes`.log
julia run_model.jl --do-not-solve --faithful2furini2016 --no-cut-position $instances 2>&1 | tee experiment_2_`date -Iminutes`.log
julia run_model.jl --do-not-solve --faithful2furini2016 $instances 2>&1 | tee experiment_3_`date -Iminutes`.log
julia run_model.jl --do-not-solve --faithful2furini2016 --round2disc $instances 2>&1 | tee experiment_17_`date -Iminutes`.log

for i in 1 2 3 4 5 6 7 8 9 10; do
	julia run_model.jl --faithful2furini2016 --no-{cut-position,redundant-cut} --seed "$i" --time-limit "$tl" $easy_for_furini 2>&1 | tee experiment_4_${i}_`date -Iminutes`.log
	julia run_model.jl --faithful2furini2016 --seed "$i" --time-limit "$tl" $easy_for_furini 2>&1 | tee experiment_5_${i}_`date -Iminutes`.log
	julia run_model.jl --faithful2furini2016 --round2disc --seed "$i" --time-limit "$tl" $easy_for_furini 2>&1 | tee experiment_6_${i}_`date -Iminutes`.log
	julia run_model.jl --seed "$i" --time-limit "$tl" $instances 2>&1 | tee experiment_7_${i}_`date -Iminutes`.log
	julia run_model.jl --round2disc --seed "$i" --time-limit "$tl" $instances 2>&1 | tee experiment_8_${i}_`date -Iminutes`.log
	julia run_model.jl --no-{cut-position,redundant-cut} --seed "$i" --time-limit "$tl" $instances 2>&1 | tee experiment_18_${i}_`date -Iminutes`.log
done

