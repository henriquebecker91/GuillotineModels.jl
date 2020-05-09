#!/bin/bash

instance_filepath="../../phd/instances/A5"
mv test_MIP_start.log test_MIP_start.log.bak
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing none --PPG2KP-MIP-start none 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing none --PPG2KP-MIP-start expected 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing none --PPG2KP-MIP-start guaranteed 2>&1 | tee -a test_MIP_start.log

./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing becker --PPG2KP-MIP-start none 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing becker --PPG2KP-MIP-start expected 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing becker --PPG2KP-MIP-start guaranteed 2>&1 | tee -a test_MIP_start.log

./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing furini --PPG2KP-MIP-start none 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing furini --PPG2KP-MIP-start expected 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-pricing furini --PPG2KP-MIP-start guaranteed 2>&1 | tee -a test_MIP_start.log

./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing none --PPG2KP-MIP-start none 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing none --PPG2KP-MIP-start expected 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing none --PPG2KP-MIP-start guaranteed 2>&1 | tee -a test_MIP_start.log

./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing becker --PPG2KP-MIP-start none 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing becker --PPG2KP-MIP-start expected 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing becker --PPG2KP-MIP-start guaranteed 2>&1 | tee -a test_MIP_start.log

./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing furini --PPG2KP-MIP-start none 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing furini --PPG2KP-MIP-start expected 2>&1 | tee -a test_MIP_start.log
./run_model.jl PPG2KP Gurobi $instance_filepath --PPG2KP-verbose --warm-jit no --PPG2KP-faithful2furini2016 --PPG2KP-pricing furini --PPG2KP-MIP-start guaranteed 2>&1 | tee -a test_MIP_start.log
