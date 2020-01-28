import GuillotineModels
import TimerOutputs
# IMPORTANT: the GuillotineModels.CommandLine.run expects a solver to be
# specified (by its julia package name, case-sensitive) in the command-line
# arguments passed to this script. Such solver needs to be imported below, or
# GuillotineModels.CommandLine.run will throw an exception. The solvers are
# not imported inside the GuillotineModels package either because then all
# solvers would become required dependencies (and users will probably have
# just a subset of such solvers installed).

# Supported solvers as of Tue 28 Jan 2020 02:06:06 AM -03:
#import CPLEX
#import Gurobi
#import Cbc # Cbc does not support --no-solver-tricks for now.

# New solvers may be added by defining a function
# empty_configured_model(::Val{:SolverName}, p_args) (which return a
# JuMP.direct_model with the specified solver already configured, see
# SolversArgs.jl file for more info), and importing the solver package
# here, as it is done with the already supported solvers.

GuillotineModels.CommandLine.run()
TimerOutputs.print_timer() # global timer used

