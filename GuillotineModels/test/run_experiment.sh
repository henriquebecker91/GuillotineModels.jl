#!/bin/bash

julia run_allsubplates_model.jl ./instances/{gcut1,gcut2,gcut3,gcut4,gcut5,gcut6,gcut7,gcut8,gcut9,gcut10,gcut11,gcut12,A5,CHL1s,CU1,CU2,CW1,CW2,CW3,CW4,CW5,CW6,CW7,CW8,CW9,CW10,CW11,okp1} | tee experiment_`date -Iminutes`.log

