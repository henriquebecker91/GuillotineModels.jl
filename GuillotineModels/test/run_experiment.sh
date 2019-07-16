#!/bin/bash

echo "ROUND 1"
julia run_allsubplates_model.jl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

echo "ROUND 2"
julia run_allsubplates_model.jl --disable-cplex-tricks ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

echo "ROUND 3"
julia run_allsubplates_model.jl --threads 2 ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

