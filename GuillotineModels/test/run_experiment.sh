#!/bin/bash

dtl="1.8E+6" # estimative of about 1 hour in the tiamats

for i in 1 2 3 4 5; do
	echo "DEFAULT, SEED $i"
	julia run_model.jl --seed "$i" --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/STS4{,s} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_${i}_`date -Iminutes`.log
done

