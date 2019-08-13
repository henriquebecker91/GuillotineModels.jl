#!/bin/bash

dtl="1.8E+6" # estimative of about 1 hour in the tiamats

for i in 1 2 3 4 5; do
	echo "WITH TRICKS SEED $i"
	julia run_model.jl --break-hvcut-symmetry --seed "$i" --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log
	echo "NO TRICKS SEED $i"
	julia run_model.jl --break-hvcut-symmetry --seed "$i" --disable-cplex-tricks --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log
done

