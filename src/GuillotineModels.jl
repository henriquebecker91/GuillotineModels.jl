module GuillotineModels

# Just include all submodules. All methods are part of some submodule and
# not of the top-level. The module names are not exported, it is expected
# the user will abbreviate the module used with a constant global variable.

include("InstanceReader.jl")
include("Utilities.jl")
include("KnapsackPlates.jl")
include("flow/Flow.jl")
include("ppg2kp/PPG2KP.jl")
include("CommandLine.jl")

end # module

