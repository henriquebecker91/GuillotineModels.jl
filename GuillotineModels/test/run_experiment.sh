#!/bin/bash

dtl="1.8E-6" # estimative of about 1 hour in the tiamats
echo "WITH TRICKS SEED 1"
julia run_allsubplates_model.jl --seed 1 --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

echo "WITH TRICKS SEED 2"
julia run_allsubplates_model.jl --seed 2 --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

echo "WITH TRICKS SEED 3"
julia run_allsubplates_model.jl --seed 3 --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

echo "NO TRICKS SEED 1"
julia run_allsubplates_model.jl --seed 1 --disable-cplex-tricks --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

echo "NO TRICKS SEED 2"
julia run_allsubplates_model.jl --seed 2 --disable-cplex-tricks --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

echo "NO TRICKS SEED 3"
julia run_allsubplates_model.jl --seed 3 --disable-cplex-tricks --det-time-limit $dtl ./instances/gcut{2,3,4,5,6,7,8,9,10,11,12} ./instances/A5 ./instances/Hchl{2,3s,4s,6s,7s} ./instances/CW{1,2,3} ./instances/CHL{1,1s,6,7} ./instances/CU{1,2} ./instances/okp{1,2,3,4,5} | tee experiment_`date -Iminutes`.log

