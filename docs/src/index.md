# GuillotineModels.jl

The main utility of this package is to enable researchers to reproduce results for the FMT formulation from [`10.1287/ijoc.2016.0710`](https://doi.org/10.1287/ijoc.2016.0710), and the enhanced formulation from [`https://arxiv.org/abs/2111.06348`](https://arxiv.org/abs/2111.06348) (final link to be updated). This can be done without knowledge of julia programming by means of the script `gmodels` which is explained in the [`README.md`](https://github.com/henriquebecker91/GuillotineModels.jl). However, the package also makes available a trove of utility functions related to the solving of the G2KP and related problems (G2MKP, G2CSP, G2OPP). Here we give a short code walkthrough from someone interested in writing Julia code that makes use of the package functions directly.

The code is divided into the following modules:

* `GuillotineModels` -- The main module. The most important types and functions pertain to this module, and are fully described at the end of this page.
* `GuillotineModels.CommandLine` -- The module used to implement the `gmodels` script, which is basically just a call to `run`. Most of the other functions in the module are only relevant to someone extending the `run` command. The only exception is `round_instance` which may be used to divide all values from a problem instance by some factor.
* `GuillotineModels.CommandLine.SolversArgs` -- The module responsible for the interface between the rest of the code and GLPK.jl, CPLEX.jl, Gurobi.jl, and Cbc.jl. It only interests someone trying to extend the set of supported solvers.
* `GuillotineModels.Utilities` -- The module aggregates the functions that are needed by one or more of other submodules but that do not are directly related to the package purpose. For example, `relax!` and `restore!` which help to change variables from binary/integer to continuous and back, or `optimize_within_time_limit!` which calls `JuMP.optimize!` but first checks for timeout and also changes the model object to respect the remaining time before timeout. The user may need to import the `SortedLinkedLW` if they want to work directly with `GuillotineModels.PPG2KP.Enumeration.gen_cuts`.
* `GuillotineModels.Utilities.Args` -- The module implements the type `Arg` widely used in `GuillotineModels.CommandLine` to represent command-line options/flags. Necessary only for understanding how to extend `GuillotineModels.CommandLine` with more options.
* `GuillotineModels.Data` -- Module responsible for all functions related to instance-reading. Useful for any user that intends to read the instances themselves, instead of relying on `GuillotineModels.CommandLine.run`. The multiple methods of function `read_from_string` are the most relevant contribution of the module. The `read_from_file` function is a convenience that just reads the whole content of a file and calls `read_from_string` over it. It also provides some limited instance writing capabilities.
* `GuillotineModels.SaveModel` -- The module provides a function `write_to_file` that saves JuMP models to some file format. The module was developed because the [`MathOptInterface.write_to_file`](https://jump.dev/MathOptInterface.jl/v0.9.19/apireference/#MathOptInterface.write_to_file) was inneficient when saving large models to MPS files. The `GuillotineModels.SaveModel.write_to_file` tries to call the annexed solver internal machinery to save the model to the desired format and, if this is not possible, it fallbacks to using the `MathOptInterface` method (which does not rely on the solver, but instead inspects the `JuMP.model` object and create the formatted output itself).
* `GuillotineModels.PPG2KP` -- The module that implements the two formulations previously mentioned. For most users, the only interesting methods are the implementations of `GuillotineModels.build_model` for the formulations mentioned.
* `GuillotineModels.PPG2KP.Heuristic` -- The module implements the 2-stage guillotine heuristic used in the pricing of [10.1287/ijoc.2016.0710](https://doi.org/10.1287/ijoc.2016.0710) (the heuristic is first presented in [10.1016/j.cor.2010.12.018](https://doi.org/10.1016/j.cor.2010.12.018)). The module can, therefore, be useful to an user interested in a simple heuristic for the G2KP for MIP-start, adding lower bound constraints, or gauging the quality of a model solution.
* `GuillotineModels.PPG2KP.Enumeration` -- The module implements the plate enumeration needed by both mainly supported formulations. The module is of little interest for the average user, but `gen_cuts` may be used to directly generate the list of all plates, cuts, and piece extractions without needing to build the model itself.
* `GuillotineModels.Flow{.Format,.Enumeration}` -- The submodule `Flow` and its submodules implement an unpublished formulation of little success. The module is kept for historical reasons. The only supported problem is G2KP and each component of an instance (lengths, widths, ...) needs to be given separatedely.

To complement the walkthrough above, we also point out the small collection of code examples in the `/examples` folder of the `GuillotineModels.jl` package.

* [`create_and_solve_PPG2KP_model.jl`](https://github.com/henriquebecker91/GuillotineModels.jl/blob/master/examples/create_and_solve_PPG2KP_model.jl) -- Create an FMT model ([`10.1287/ijoc.2016.0710`](https://doi.org/10.1287/ijoc.2016.0710)) for a small instance and solve it, without using `GuillotineModels.CommandLine.run`.
* [`instance_converter.jl`](https://github.com/henriquebecker91/GuillotineModels.jl/blob/master/examples/instance_converter.jl) -- Example on how to extend `GuillotineModels.Data.write_to_file`. Reads a `Classic_G2KP` instance file and save it in the `CPG_SLOPP` format.
* [`run_discretization.jl`](https://github.com/henriquebecker91/GuillotineModels.jl/blob/master/examples/run_discretization.jl) -- Calls `GuillotineModels.PPG2KP.Enumeration.gen_cuts` directly and print some info about the descretization of a small instance.
* [`run_heuristic.jl`](https://github.com/henriquebecker91/GuillotineModels.jl/blob/master/examples/run_heuristic.jl) -- Takes an instance in the `Classic_G2KP` format and run both the unoptimized and the optimized versions of the Guillotine 2-stage heuristic from `GuillotineModels.PPG2KP.Heuristic` over it. Shows the solution and the timings.

For a better understanding on how to deal with already built models we suggest reviewing the [`JuMP documentation for version 0.21.5`](https://jump.dev/JuMP.jl/v0.21.5/).

The code heavily relies on parametric types, both for templated type and method signatures. Most of the cases, the type parameters are named `D`, `S` , `P`. The three are expected to be concrete subtypes of `Integer`. The `D` represents the type used to store piece *demand* (should be large enough to store the sum of all pieces demand). The `S` represents the type used to store the pieces *size* (should be enough to store the sum of all pieces length or width). The `P` represents the type used to store *profit* (or area, should be large enough to store the original plate area or the sum of all piece profits, whichever is the largest).

```@meta
CurrentModule = GuillotineModels
```

```@autodocs
Modules = [GuillotineModels]
```

