#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "usage: ./script.sh log_with_prints"
    exit 1
fi

logfile="$1"

csv=""
header="instfname;seed;n;n_;p_;time_to_enumerate_plates;time_to_solver_build;time_to_build_model;time_to_solve_model;num_vars;num_constrs;num_plates;num_picuts;num_cuts_made;obj_value;obj_bound;stop_code"
IFS=';' read -ra aheader <<< "$header"
for f in "${aheader[@]}"; do
	column=`grep "^$f = " "$logfile" | cut -f 3- -d\  `
	if [ -z "$csv" ]; then
		csv="$column"
	else
		csv=`paste -d \; <(echo "$csv") <(echo "$column")`
	fi
done
echo "$header"
echo "$csv"

