# Check if the necessary packages are installed and define the project
# to make them available to importation by "using".
#using Pkg
#Pkg.activate("..")
#Pkg.instantiate()

import GuillotineModels
import TimeOutputs

#TimerOutputs.enable_debug_timings(GuillotineModels)
GuillotineModels.RunModel.run_batch()
TimeOutputs.print_timer() # global timer used

