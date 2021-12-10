# GuillotineModels.jl

Mathematical models for 2D guillotine cutting problems. Part of Henrique Becker's PhD (co-advised by Olinto Araujo, advised by Luciana S. Buriol).

The best maintened formulations are:
1. FMT from https://doi.org/10.1287/ijoc.2016.0710, with and without pricing.
2. An enhanced formulation based on FMT. For now only the preprint is available: https://arxiv.org/abs/2111.06348, but the link will be updated when the paper is accepted.

The code documentation can be found in https://henriquebecker91.github.io/GuillotineModels.jl/stable/

## Installation and usage tutorial

As the code may be of interest of academics that do not have familiarity with the Julia programming language, the next steps also describe the usual steps to get ay Julia package working.

### Julia installation

You can skip these steps if you have Julia (1.4.2 or superior) installed.

The recommended version of Julia is 1.4.2, which can be found in https://julialang.org/downloads/oldreleases/. Given how semantic versioning works (https://semver.org/), the code should work with any 1.y.z versions in which `y >= 4`. Julia 1.4.2 is just the version in which the code was more extensively tested.

The download page provides installers for Windows and Mac. Note, however, the code was only tested in Linux. In Linux, it is necessary to unpack the downloaded tarball, add the folder `julia-1.4.2/bin` to the PATH and, optionally, automate the change of `PATH` to happen when the shell of your preference starts. If you use Bash, then you can do:

```
wget https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz
tar -xf julia-1.4.2-linux-x86_64.tar.gz
export PATH="`pwd`/julia-1.4.2/bin:$PATH"
echo -e "\nexport PATH=\"`pwd`/julia-1.4.2/bin:\$PATH\"" >> $HOME/.bashrc
```

The commands assume a 64-bits x86 glibc Linux and may need adaptation if the assumption is wrong. For a different version, or if the link is broken, check the Julia Language downloads page (https://julialang.org/downloads/).

If the folder changes place after the installation, the PATH must changed to point to the new location of the `julia-1.4.2/bin` folder.

If everything is correct, then calling `julia` in a shell should show you a banner and a `julia>` prompt, in which `julia` code can be executed.

### Package installation

The julia code to install the package follows:

```
import Pkg
Pkg.add(Pkg.PackageSpec(url = "https://github.com/henriquebecker91/GuillotineModels.jl"))
```

In Linux, the package is installed at `~/.julia/packages/GuillotineModels/`. This folder will have one or more subdirectories (one for each distinct installed version of the package). For now it should have just one, but if distinguishing between multiple installed versions becomes necessary, then executing the following command from inside the folder should help:

```
$ grep version */Project.toml # Gives the version inside each folder.
RtyzT/Project.toml:version = "0.4.0"
Sft3a/Project.toml:version = "0.4.2"
wEW7Q/Project.toml:version = "0.6.5"
xtvPM/Project.toml:version = "0.7.0"
# ...
```

### Command-line usage

The best way to use the package to solve instances with minimal julia knowledge is to call the `gmodels` script. The script is inside the `scripts` folder of the package. The previous section has details on how to find the package folder. The first time the script runs it will install any other dependencies, including GLPK. GLPK is the default solver because it is automatically installed by the Julia package manager when the package `GLPK.jl` is installed. In Linux/Bash, you may save the list of parameters with the following command:

```
./gmodels --help | tee help.txt
```

For convenience, the same list is annexed at the end of this README.

If you did not add the folder with the `julia` binary to the PATH variable, then it may be necessary to call the command the following way:

```
/path/to/julia --project=. ./gmodels --help | tee help.txt
```

If both the folder with `julia` and the folder with `gmodels` are added to the PATH, then `gmodels` can be called from anywhere. Executing the following Bash commands inside the `scripts` folder will add it to the PATH for the current and new Bash sessions:

```
export PATH="`pwd`:$PATH"
echo -e "\nexport PATH=\"`pwd`:\$PATH\"" >> $HOME/.bashrc
```

To solve a small G2KP (the knapsack variant) instance, create a plain text file with the following content:

```
250 250
10
167 184 30728 1
114 118 13452 1
167 152 25384 1
83 140 11620 1
70 86 6020 1
143 166 23738 1
120 160 19200 1
66 148 9768 1
87 141 12267 1
69 165 11385 1
```

If you named the file `instance.txt` and it is in the same folder as gmodels, then you may solve it with the FMT model (https://doi.org/10.1287/ijoc.2016.0710) by calling:

```
./gmodels G2KP Classic_G2KP PPG2KP GLPK ./instance.txt --PPG2KP-faithful2furini2016 | tee -a output.txt
```

To solve an instance with the enhanced formulation developed by the authors of this repository, you may call:

```
./gmodels G2KP Classic_G2KP PPG2KP GLPK ./instance.txt --PPG2KP-round2disc | tee -a output.txt
```

By default, no pricing procedure is applied. Adding `--PPG2KP-pricing furini` to the command should enable the pricing descibed in (https://doi.org/10.1287/ijoc.2016.0710). It works on both the FMT formulation and the enhanced one.

### Reading the output

The solving process output quite a lot of information by default. All time durations are in seconds (unless acompannied by a distinct suffix). Also, given how Julia works, the timings of each section/function from `gmodels` output include compilation, to avoid this check the parameter `--warm-jit` explained by `--help`. The output of the first solving command of the last section is explained below:

```
 Activating environment at `~/Aulas/doutorado/GuillotineModels.jl/scripts/Project.toml`
 ...
```

Any messages before `args = ...` are related to the code that checks if every package needed installed. After the first run, in which the packages are installed, this should output only the line above.

```
args = ["G2KP", "Classic_G2KP", "PPG2KP", "GLPK", "./instance.txt"]
```

The arguments passed to the command in the form of a julia `Vector{String}` literal. Easy to copy to a julia code.

```
p_args = Dict{String,Any}("PPG2KP-verbose" => false,"save-model" => "","PPG2KP-no-redundant-cut" => false,"relax2lp" => false,"PPG2KP-pricing-beta" => 0.25,"PPG2KP-no-furini-symmbreak" => false,"PPG2KP-heuristic-seed" => 1,"PPG2KP-quiet" => false,"use-native-save-model-if-possible" => false,"PPG2KP-ignore-2th-dim" => false,"GLPK-raw-parameters" => "Pair{String, Any}[]","PPG2KP-building-time-limit" => 3.1536e7,"instance_path" => "./instance.txt","PPG2KP-hybridize-with-restricted" => false,"PPG2KP-ignore-d" => false,"PPG2KP-aggressive-hybridization" => false,"format" => "Classic_G2KP","no-csv-output" => false,"problem" => "G2KP","warm-jit" => "no","PPG2KP-do-not-purge-unreachable" => false,"model" => "PPG2KP","PPG2KP-pricing" => "none","GLPK-no-output" => false,"PPG2KP-no-cut-position" => false,"do-not-solve" => false,"solver" => "GLPK","round-up" => 1.0,"PPG2KP-mirror-plates" => false,"PPG2KP-pricing-alpha" => 0.2,"PPG2KP-use-c25" => false,"GLPK-time-limit" => 2097152,"round-nearest" => 1.0,"round-down" => 1.0,"PPG2KP-faithful2furini2016" => false,"PPG2KP-Gurobi-LP-method-inside-furini-pricing" => -2,"PPG2KP-MIP-start" => "expected","PPG2KP-minimal-subset-in-iterative-pricing" => false,"PPG2KP-no-sort-in-iterative-pricing" => false,"generic-time-limit" => 3.1536e7,"PPG2KP-allow-rotation" => false,"PPG2KP-round2disc" => false)
```

The julia dictionary literal containing the value used for every parameter, even the ones not explicitly given. While the default values can be queried with `--help` this can be useful because a parameter may change default between versions (so this enhances reproducibility).

```
date_now = "2021-12-03T15:42:07"
instance_path = "./instance.txt"
read_instance_time = 0.028048992156982422
```

The date which the script started, the instance path, and the time (in seconds) used to read the instance(s). Some accepted formats allow for more than a single instance inside a single file. These information are not instance-specific (i.e., related to the whole file). After them start the instance-specific information.

```
INSTANCE_START_MARKER_1
L = 250
W = 250
n_ = 10
n = 10
min_l = 66
min_w = 86
p_ = 3.787878787878788
GuillotineModels.Data.G2KP{Int64,Int64,Int64}(250, 250, [167, 114, 167, 83, 70, 143, 120, 66, 87, 69], [184, 118, 152, 140, 86, 166, 160, 148, 141, 165], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [30728, 13452, 25384, 11620, 6020, 23738, 19200, 9768, 12267, 11385])
```

Each distinct instance inside the same file is enclosed in `INSTANCE_START_MARKER_X` and `INSTANCE_CLOSE_MARKER_X` where `X` is the number of the instance inside the file. This happens even if the file has a single instance inside. The above output prints basic properties of the instance with their default names. The only ones which merit an explanation are: `n_` -- sum of the piece demand (distinct from `n` which is the number of piece types), and `p_` -- the largest ratio between a dimension of a piece and the corresponding dimension in the original plate.

```
build_stop_reason = BUILT_MODEL
build_stop_code = 0
num_vars = 1052
num_constrs = 3557
build_time = 6.362917900085449
```

This block brings information associated with the model generation. The `build_stop_reason = BUILT_MODEL` means the model construction finished with the built model; it may also return `FOUND_OPTIMUM` when pricing is used, because the pricing procedure may find and prove optimality before finishing to build the model.

```
MARK_FINAL_GENERIC_SOLVE
*     0: obj =  -0.000000000e+00 inf =   0.000e+00 (669)
*   479: obj =   5.073950000e+04 inf =   1.110e-16 (0)
+   479: mip =     not found yet <=              +inf        (1; 0)
Solution found by heuristic: 48368
+   530: mip =   4.836800000e+04 <=     tree is empty   0.0% (0; 25)
finished_model_solve = 0.109340028
stop_reason = MathOptInterface.OPTIMAL
stop_code = 1
build_and_solve_time = 6.5119640827178955
obj_value = 48368.0
obj_bound = 48368.0
```

The `MARK_FINAL_GENERIC_SOLVE` marks the start of the solver output for the final model, this is, the model after any pricing procedures. It is followed by some information on the solving process: `stop_reason` and `stop_code` -- the reason for the solver stopping, the list of possible values is available at https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.TerminationStatusCode; `build_and_solve_time` -- the time to build the final model (including pricing) and to solve it; `obj_value` -- the objective function value of the incumbent solution; `obj_bound` -- the best bound on the objective function value found during the solving process.

```
solution = GuillotineModels.CutPattern{Int64,Int64}(250, 250, 0, false, GuillotineModels.CutPattern{Int64,Int64}[GuillotineModels.CutPattern{Int64,Int64}(83, 250, 0, true, GuillotineModels.CutPattern{Int64,Int64}[GuillotineModels.CutPattern{Int64,Int64}(83, 86, 0, false, GuillotineModels.CutPattern{Int64,Int64}[GuillotineModels.CutPattern{Int64,Int64}(70, 86, 5, false, GuillotineModels.CutPattern{Int64,Int64}[])]), GuillotineModels.CutPattern{Int64,Int64}(83, 164, 0, false, GuillotineModels.CutPattern{Int64,Int64}[GuillotineModels.CutPattern{Int64,Int64}(83, 140, 4, false, GuillotineModels.CutPattern{Int64,Int64}[])])]), GuillotineModels.CutPattern{Int64,Int64}(167, 250, 0, false, GuillotineModels.CutPattern{Int64,Int64}[GuillotineModels.CutPattern{Int64,Int64}(167, 184, 1, false, GuillotineModels.CutPattern{Int64,Int64}[])])])
```

The `solution = ...` block has the literal `CutPattern` object that may be used to input the solution back to a julia session (given package `GuillotineModels` is imported). This solution is the incumbent solution. Besides this representation, many alternative representations of this same solution are printed by the code.

```
PRETTY_STR_SOLUTION_BEGIN
P250x250{
  P83x250[
    P83x86{ 5p70x86 }
    P83x164{ 4p83x140 }
  ]
  P167x250{ 1p167x184 }
}
PRETTY_STR_SOLUTION_END
```

The representation above may be obtained by passing the `solution = ...` object to the [`GuillotineModels.to_pretty_str`](https://henriquebecker91.github.io/GuillotineModels.jl/stable/#GuillotineModels.to_pretty_str-Union{Tuple{CutPattern{D,S}},%20Tuple{S},%20Tuple{D}}%20where%20S%20where%20D) function. The explanation of the representation can be found in the documentation of the aforementioned function.

```
TIKZ_PRETTY_STR_SOLUTION_BEGIN
\begin{tikzpicture}
\draw[dashed, thick, black] (0, 0) rectangle (70, 86);
\node[font=\LARGE] at (35.0, 43.0) {5};
%\node[font=\LARGE] at (35.0, 43.0) {70x86};
\draw[dashed, thick, black] (0, 0) rectangle (83, 86);
\draw[dashed, thick, black] (0, 86) rectangle (83, 226);
\node[font=\LARGE] at (41.5, 156.0) {4};
%\node[font=\LARGE] at (41.5, 156.0) {83x140};
\draw[dashed, thick, black] (0, 86) rectangle (83, 250);
\draw[dashed, thick, black] (0, 0) rectangle (83, 250);
\draw[dashed, thick, black] (83, 0) rectangle (250, 184);
\node[font=\LARGE] at (166.5, 92.0) {1};
%\node[font=\LARGE] at (166.5, 92.0) {167x184};
\draw[dashed, thick, black] (83, 0) rectangle (250, 250);
\draw[dashed, thick, black] (0, 0) rectangle (250, 250);
\end{tikzpicture}
TIKZ_PRETTY_STR_SOLUTION_END
```

The representation above can be used to obtain a Tikz diagram in a LaTex document. The authors suggest the `qtikz` application for easily checking the diagram without having to embed it to a complete LaTex document. The function used to obtain this representation from the `CutPattern` object is [`GuillotineModels.to_tikz_picture`](https://henriquebecker91.github.io/GuillotineModels.jl/stable/#GuillotineModels.to_tikz_picture-Union{Tuple{CutPattern{D,S}},%20Tuple{S},%20Tuple{D}}%20where%20S%20where%20D).

```
SIMPLIFIED_PRETTY_STR_SOLUTION_BEGIN
P250x250{
  P83x250[ 5p70x86 4p83x140 ]
  1p167x184
}
SIMPLIFIED_PRETTY_STR_SOLUTION_END
TIKZ_SIMPLIFIED_PRETTY_STR_SOLUTION_BEGIN
\begin{tikzpicture}
\draw[dashed, thick, black] (0, 0) rectangle (70, 86);
\node[font=\LARGE] at (35.0, 43.0) {5};
%\node[font=\LARGE] at (35.0, 43.0) {70x86};
\draw[dashed, thick, black] (0, 86) rectangle (83, 226);
\node[font=\LARGE] at (41.5, 156.0) {4};
%\node[font=\LARGE] at (41.5, 156.0) {83x140};
\draw[dashed, thick, black] (0, 0) rectangle (83, 250);
\draw[dashed, thick, black] (83, 0) rectangle (250, 184);
\node[font=\LARGE] at (166.5, 92.0) {1};
%\node[font=\LARGE] at (166.5, 92.0) {167x184};
\draw[dashed, thick, black] (0, 0) rectangle (250, 250);
\end{tikzpicture}
TIKZ_SIMPLIFIED_PRETTY_STR_SOLUTION_END
```

The block above shows both representations already described, the only distinction is that [`GuillotineModels.simplify!`](https://henriquebecker91.github.io/GuillotineModels.jl/stable/#GuillotineModels.simplify!-Union{Tuple{CutPattern{D,S}},%20Tuple{S},%20Tuple{D}}%20where%20S%20where%20D) was applied to the `CutPattern` object before printing. The function rearranges the `CutPattern` object to remove unnecessary cuts sometimes present in the model. Most formulations allow further cutting waste plates without necessity.

```
solution_value = 48368
solution_print_time = 2.412231922149658
total_instance_time = 12.945231914520264
INSTANCE_CLOSE_MARKER_1
```

Above there is the last block of information about the specific instance. The `solution_value` is computed based solely on the `CutPattern` object, and should match the solver `obj_value`.

```
total_file_time = 13.264692852
```

As mentioned before, some formats allow a single file to contain multiple instances, the above is the total aggregated time to solve all instances.

```
 ──────────────────────────────────────────────────────────────────
                                                     Time          
                                             ──────────────────────
              Tot / % measured:                   44.3s / 38.1%    

 Section                             ncalls     time   %tot     avg
 ──────────────────────────────────────────────────────────────────
 run                                      1    16.9s   100%   16.9s
   read_build_solve_and_print             1    13.0s  76.9%   13.0s
     build_model                          1    3.20s  19.0%   3.20s
       _build_base_model!                 1    983ms  5.82%   983ms
       _gen_cuts_wo                       1    819ms  4.85%   819ms
       _handle_unreachable!               1    516ms  3.06%   516ms
         _purge_unreachable!              1   3.72ms  0.02%  3.72ms
         _reachable                       1   67.8μs  0.00%  67.8μs
           _reachable_plate_types         1   63.0μs  0.00%  63.0μs
           _reachable_carves              2   2.52μs  0.00%  1.26μs
         _print_unreachable_plates        1    419ns  0.00%   419ns
       VarInvIndexes                      1    128μs  0.00%   128μs
     get_cut_pattern                      1    1.01s  5.98%   1.01s
     stringfy_solutions                   1    554ms  3.28%   554ms
     finished_model_solve                 1    109ms  0.65%   109ms
     empty_configured_model               1   75.0ms  0.44%  75.0ms
     round_instance                       1   27.1ms  0.16%  27.1ms
     read_from_file                       1   15.5ms  0.09%  15.5ms
       read_from_string                   1   15.4ms  0.09%  15.4ms
     print_solutions                      1   28.1μs  0.00%  28.1μs
   parse_args                             1    3.37s  19.9%   3.37s
   warm_jit                               1   3.61μs  0.00%  3.61μs
 ──────────────────────────────────────────────────────────────────
```

Finally, the table above shows how much of the total time was spent on some key procedures. Again, if `--warm-jit` was not changed from the default value, then each function includes the compilation time of the functions it calls. In the example above, the compilation time explains most of the disparity between the procedure total times and the sum of their parts. If the input file contained multiple instances then this table aggregates the information for all solves.

### Installing other solver for use

Besides GPLK, the code aims to support CPLEX, Gurobi, and Cbc. It is also possible to extend the tool to support other solvers without having to change the package but this requires some julia knowledge (check [`src/SolverArgs.jl`](https://github.com/henriquebecker91/GuillotineModels.jl/blob/master/src/SolversArgs.jl) for more details).

There is a julia package for each of these solvers. They take care of the communication with the solver. The most recent version of these packages for which the code was tested are: GLPK.jl -- v0.14.14, CPLEX.jl -- v0.7.8, Gurobi.j -- v0.9.14, and Cbc.jl -- v0.8.1. The version of the julia package is not the same as the version of the solver itself. Both GLPK.jl (already configured) and Cbc.jl automatically install version v4.64.0 (GLPK) and v2.10.3 (Cbc); but they may be used with custom installations if desired (check out their GitHub pages for instructions: https://github.com/jump-dev/GLPK.jl/tree/v0.14.14 and https://github.com/jump-dev/Cbc.jl/tree/v0.8.1). The CPLEX.jl and Gurobi.jl packages need you to independently install the solver. The supported solver versions by the tested package versions are 9.0--9.1 (Gurobi) and 12.10--20.1 (CPLEX).

The procedure for changing the `gmodels` script to support the solvers is described below.

1. (Gurobi and CPLEX only) **Obtain a legal copy of the commercial solver.** The GitHub page of the julia package (of the correct version) corresponding to the desired solver (https://github.com/jump-dev/Gurobi.jl/tree/v0.9.14 or https://github.com/jump-dev/CPLEX.jl/tree/v0.7.8) has, in its README, links the official website of the solvers, and also give some advice on how to legally obtain a copy. Take note of the path in which the solver was installed.
2. **Open the REPL for installing the julia package.** Execute `julia --project=.` inside the `scripts` folder of the `GuillotineModels` package. See sections *Package installation* and *Command-line usage* of this README for more info on how to find this folder.
3. (Gurobi and CPLEX only) **Set the path of the solver.** As pointed out in the already linked GitHub page of the Gurobi.jl and CPLEX.jl packages, it is necessary to set either `ENV["CPLEX_STUDIO_BINARIES"] = "/path/to/os_arch/folder/inside/bin"` or `ENV["GUROBI_HOME"] = "/path/to/os_arch/folder"` in the Julia REPL before installing the package. Check the GitHub page instructions if you are not sure exactly which folder the path should be pointing at.
4. **Install the Julia package corresponding to the solver.** In the same Julia REPL session execute: `import Pkg; Pkg.add("SOLVER_NAME"); Pkg.build("SOLVER_NAME")`. If the build process fails then probably the path set to `ENV` in the last step is not correct, only `Pkg.build("SOLVER_NAME")` will need to be re-run after fixing `ENV`. The command `Pkg.status()` may be used to check the version of the installed package.
5. **Edit the `gmodels` script.** Finally, the `gmodels` script needs to be edited with a plain text editor (some examples are `nano` or `gedit` on Linux, `TextEdit` on Mac, and `Notepad++` on Windows). The line `#import SOLVER_NAME` needs to have the `#` removed.

### The instance formats

The code initially supported only the G2KP (Guillotine 2D Knapsack Problem) and the simplest instance format used for this problem, which we call `Classic_G2KP` and is described in the documentation [`here`](https://henriquebecker91.github.io/GuillotineModels.jl/stable/Data/#GuillotineModels.Data.read_from_string-Union{Tuple{P},%20Tuple{S},%20Tuple{D},%20Tuple{GuillotineModels.Data.Classic_G2KP{D,S,P},AbstractString}}%20where%20P%20where%20S%20where%20D).

The code now supports many formats from the 2DCPackGen ([10.1016/j.ejor.2014.02.059](https://doi.org/10.1016/j.ejor.2014.02.059), see [https://sites.google.com/gcloud.fe.up.pt/cutting-and-packing-tools/2dcpackgen](https://sites.google.com/gcloud.fe.up.pt/cutting-and-packing-tools/2dcpackgen) for the code). Each of these formats starts with a self-explanatory header, which is present in the documentation of the specific [`read_from_string`](https://henriquebecker91.github.io/GuillotineModels.jl/stable/Data/#GuillotineModels.Data.read_from_string-Union{Tuple{P}, Tuple{S}, Tuple{D}, Tuple{GuillotineModels.Data.CPG_Format{D,S,P},AbstractString}} where P where S where D) method.

Besides the G2KP, the code now also support the G2MKP (Guillotine 2D Multiple Knapsack Problem), G2OPP (Guillotine 2D Orthogonal Packing Problem, i.e., the decision problem), and G2CSP (Guillotine 2D Cutting Stock Problem). Only the homogeneous variant (i.e., all original plates of the same size) of G2MKP and of G2CSP are supported. The supported problem-format pairs are: `G2KP`-`Classic_G2KP`, `G2KP`-`CPG_SLOPP`, `G2KP`-`Simple_CPG_SLOPP`, `G2KP`-`CPG_SSSCSP`, `G2CSP`-`CPG_SSSCSP`, `G2CSP`-`Classic_G2KP`, `G2OPP`-`CPG_SSSCSP`, `G2MKP`-`CPG_MHLOPPW`.

### Output of `--help`

```
usage: gmodels [--do-not-solve]
               [--generic-time-limit GENERIC-TIME-LIMIT]
               [--save-model SAVE-MODEL]
               [--use-native-save-model-if-possible]
               [--warm-jit WARM-JIT] [--no-csv-output] [--relax2lp]
               [--round-nearest ROUND-NEAREST] [--round-up ROUND-UP]
               [--round-down ROUND-DOWN]
               [--GLPK-time-limit GLPK-TIME-LIMIT] [--GLPK-no-output]
               [--GLPK-raw-parameters GLPK-RAW-PARAMETERS]
               [--PPG2KP-minimal-subset-in-iterative-pricing]
               [--PPG2KP-no-sort-in-iterative-pricing]
               [--PPG2KP-building-time-limit PPG2KP-BUILDING-TIME-LIMIT]
                              [--PPG2KP-Gurobi-LP-method-inside-furini-pricing PPG2KP-GUROBI-LP-METHOD-INSIDE-FURINI-PRICING]
               [--PPG2KP-pricing-alpha PPG2KP-PRICING-ALPHA]
               [--PPG2KP-pricing-beta PPG2KP-PRICING-BETA]
               [--PPG2KP-round2disc]
               [--PPG2KP-hybridize-with-restricted]
               [--PPG2KP-faithful2furini2016]
               [--PPG2KP-no-redundant-cut] [--PPG2KP-no-cut-position]
               [--PPG2KP-no-furini-symmbreak]
               [--PPG2KP-pricing PPG2KP-PRICING] [--PPG2KP-ignore-d]
               [--PPG2KP-use-c25] [--PPG2KP-ignore-2th-dim]
               [--PPG2KP-heuristic-seed PPG2KP-HEURISTIC-SEED]
               [--PPG2KP-verbose] [--PPG2KP-quiet]
               [--PPG2KP-do-not-purge-unreachable]
               [--PPG2KP-MIP-start PPG2KP-MIP-START]
               [--PPG2KP-allow-rotation] [--PPG2KP-mirror-plates]
               [--PPG2KP-aggressive-hybridization] [-h]
               [problem] [format] [model] [solver] [instance_path]

optional arguments:
  -h, --help            show this help message and exit

Core Parameters:
  problem               The type of problem to be solved (case
                        sensitive, ex.: G2KP, G2CSP, G2ODP, G2MKP).
                        Required.
  format                The format of the instance(s) in the file
                        (case sensitive, ex.: Classic_G2KP, CPG_SLOPP,
                        CPG_SSSCSP, CPG_ODPW, CPG_MHLOPPW). Required.
  model                 Model or solution procedure to be used (case
                        sensitive, ex.: Flow, KnapsackPlates, PPG2KP).
                        Required.
  solver                Solver to be used if necessary
                        (case-sensitive, ex.: NoSolver, Cbc, CPLEX,
                        Gurobi). Required, even if --do-not-solve is
                        specified. NoSolver use `JuMP.Model()` which
                        allows building and saving the model, but not
                        solving it.
  instance_path         The path to the instance to be solved.

Generic Options:
  --do-not-solve        The model is built but not solved. A solver
                        has yet to be specified. Note that just
                        building a model may depend on using a solver
                        over subproblems. Such uses of the solver are
                        not disabled by this flag.
  --generic-time-limit GENERIC-TIME-LIMIT
                        Defines a time limit (in seconds) to be
                        observed in the context of the model-agnostic
                        process of reading, building, solving, and
                        printing. Each model has to define its own
                        flag to control the time inside the model
                        building process (to be set independently, it
                        is does not interact with this flag).
                        `JuMP.set_time_limit_sec` is called over the
                        model before starting to solve it, with the
                        remaining time after reading instance and
                        building the model (not the total time). If
                        `--warm-jit` defines that there will be a mock
                        run to warm the jit, the time limit applies to
                        this mock run (i.e., if the mock run breaks
                        the limit it an exception is raised) but the
                        timer is reset after the mock, so the mock
                        time is not counted for the 'real' run time
                        limit. A `GuillotineModels.TimeoutError` is
                        raised if the time limit is violated. (type:
                        Float64, default: 3.1536e7)
  --save-model SAVE-MODEL
                        Save the model of the passed instance to the
                        given filename. The format used depends on the
                        extension of the filename, allowed formats are
                        described by the enumeration
                        `MathOptInterface.FileFormats.FileFormat` from
                        package `MathOptInterface`). The filename also
                        allows using $<parameter_name> patterns inside
                        the name, these will be replaced by the value
                        of the parameter (if it was not passed, the
                        default value will be used). For convenience,
                        a pseudo-parameter instance_file is also
                        available (it is the same as
                        `basename(instance_path)`, it includes any
                        file extensions), and `Bool` parameters (i.e.,
                        flags with no argument) are replaced by 0 or
                        1. Putting an @ in front of $<some_string>
                        will disable the substitution (and remove the
                        @ from the final string). There is no way to
                        express an @ followed by a substitution that
                        actually occurs (this is a limitation of the
                        code). (default: "")
  --use-native-save-model-if-possible
                        Only used if a filename is passed to
                        '--save-model'. If
                        `--use-native-save-model-if-possible` is
                        passed, then
                        `GuillotineModels.SaveModel.write_to_file` is
                        used instead of `JuMP.write_to_file`. Check
                        the `SaveModel` module documentation for more
                        info. For at least Gurobi and CPLEX,
                        `GuillotineModels.SaveModel.write_to_file`
                        will use the raw solver model-saving utilities
                        instead of the `JuMP` one, but if the solver
                        is not recognized then it will default to the
                        JuMP one.
  --warm-jit WARM-JIT   There are three valid options: 'no',
                        'with-toy', and 'yes'. If 'yes', the instance
                        is solved two times, but the first time the
                        arguments are changed to disable output (i.e.,
                        '--no-csv-output' and '--SOLVER-no-output' are
                        passed, and if supported, '--MODEL-quiet' is
                        passed too; in the GuillotineModels.TIMER the
                        timings of this first run are under
                        'warm_jit'). If 'with-toy', it behaves
                        similarly to 'yes' but instead of using the
                        instance itself it uses a very small hardcoded
                        instance (this helps a lot, but many
                        procedures only called when the model has
                        specific properties are not called, and
                        therefore, not compiled). If 'no', the code is
                        just run a single time. (default: "no")
  --no-csv-output       Disable some extra output that the author of
                        this package uses to assemble CSVs and/or
                        debug. Also disables --save-model.
  --relax2lp            Integer and binary variables become
                        continuous.
  --round-nearest ROUND-NEAREST
                        Multiplies the instances size-related fields
                        by the passed factor and round them to nearest
                        (the model answer becomes a GUESS, not a valid
                        primal heuristic, nor a valid bound). (type:
                        Float64, default: 1.0)
  --round-up ROUND-UP   Multiplies the instances size-related fields
                        by the passed factor and round them up (the
                        model becomes a PRIMAL HEURISTIC). (type:
                        Float64, default: 1.0)
  --round-down ROUND-DOWN
                        Multiplies the instance size-related fields by
                        the passed factor and round them down (the
                        model becomes an OPTIMISTIC GUESS, A VALID
                        BOUND). (type: Float64, default: 1.0)

GLPK-specific Options:
  --GLPK-time-limit GLPK-TIME-LIMIT
                        BROKEN, DO NOT USE, ALWAYS OVERWRITTEN BY
                        `--generic-time-limit`. Set GLPK parameter:
                        tm_lim. The original parameter is in
                        milliseconds, but to keep it similar to the
                        other solvers this option takes seconds. To
                        set this parameter with milliseconds precision
                        use the --raw-parameter option. The default is
                        a little over 24 days because GLPK uses a
                        Int32 for milliseconds and 24 days is close to
                        the maximum time-limit supported. (type:
                        Int64, default: 2097152)
  --GLPK-no-output      Set msg_lev to GLPK.MSG_OFF. Disables the
                        solver output.
  --GLPK-raw-parameters GLPK-RAW-PARAMETERS
                        A string of Julia code to be evaluated to GLPK
                        parameters. For example: 'Pair{String,
                        Any}["msg_lev" => GLPK.MSG_OFF]' will have the
                        same effect as --GLPK-no-output. Obviously,
                        this is a security breach. If you find the
                        complete list of GLPK parameters please send
                        it to me (henriquebecker91@gmail.com).
                        (default: "Pair{String, Any}[]")

PPG2KP-specific Options:
  --PPG2KP-minimal-subset-in-iterative-pricing
                        Only used if the value of the 'pricing' option
                        is 'furini'. The explanation of the iterative
                        pricing algorithm in 10.1287/ijoc.2016.0710
                        implies all rows are present since start, and
                        that cuts (`x` variables) that cannot assume
                        non-zero values yet (because the plates they
                        are cutting are not made available by other
                        cuts yet) are also present since start. This
                        flag enables an alternative behavior, in which
                        the only plate constraints in the initial
                        model are either: obtainable by restricted
                        cuts (i.e., present in a restricted version of
                        the model); or have a sell/extraction (`y`
                        variable in original model or `e` variable in
                        the revised model) selling that plate (or
                        extracting a piece from it). All
                        selling/extraction variables are present, but
                        the cut variables are a smaller subset (only
                        the cuts over 'reachable' plates, no cuts over
                        plates that cannot be yet generated by the
                        model). Consequently, both 'cut
                        variables'/columns and 'plate
                        constraints'/rows are added
                        dinamically/iteratively to the model.
  --PPG2KP-no-sort-in-iterative-pricing
                        Only used if the value of the 'pricing' option
                        is 'furini'. The explanation of the pricing
                        algorithm in 10.1287/ijoc.2016.0710 does not
                        make very clear if the variables with positive
                        reduced profit added each iteration are the
                        ones with highest reduced profit or an
                        arbitrary subset. The default is to sort
                        because: the extra effort is minimal compared
                        to solving the LP; the time solving the LP can
                        be drastically reduced by the adding the
                        variables with highest reduced profit first.
  --PPG2KP-building-time-limit PPG2KP-BUILDING-TIME-LIMIT
                        Defines a time limit (in seconds) to be
                        observed in the context of the model building
                        proccess. If the solver is called during the
                        building procedure then
                        `JuMP.set_time_limit_sec` is called over the
                        model before starting to solve it, with the
                        remaining time (not the total time). A
                        `GuillotineModels.TimeoutError` is raised if
                        the time limit is violated. (type: Float64,
                        default: 3.1536e7)
  --PPG2KP-Gurobi-LP-method-inside-furini-pricing PPG2KP-GUROBI-LP-METHOD-INSIDE-FURINI-PRICING
                        Allows to switch the algorithm used to solve
                        continuous models and root node relaxations
                        inside the Furini's iterative pricing when
                        Gurobi is the solver used. Will have no effect
                        but print warnings if either Furini's pricing
                        is not called or Gurobi is not the solver in
                        use. The default value of -2 will not touch
                        Gurobi's Method parameter, values -1 to 5 will
                        set Gurobi's Method to the corresponding value
                        during the subprocedure mentioned, and restore
                        the previous value after it. (type: Int64,
                        default: -2)
  --PPG2KP-pricing-alpha PPG2KP-PRICING-ALPHA
                        Used to compute the number of variables added
                        in each iteration of the iterative pricing
                        (only used if furini pricing is selected, see
                        pricing flag). Must be above zero and at most
                        one. Explained in 10.1287/ijoc.2016.0710, p.
                        13 (747) (last paragraph before section 4.3).
                        Default value is the one used in the
                        experiments of the original paper. (type:
                        Float64, default: 0.2)
  --PPG2KP-pricing-beta PPG2KP-PRICING-BETA
                        Used to compute the number of variables addded
                        in each iteration of the iterative pricing
                        (only used if furini pricing is selected, see
                        pricing flag). Must be non-negative, and makes
                        most sense to be at most one. Explained in
                        10.1287/ijoc.2016.0710, p. 13 (747) (last
                        paragraph before section 4.3). Default value
                        is the one used in the experiments of the
                        original paper. (type: Float64, default: 0.25)
  --PPG2KP-round2disc   Round the second child size to a discretized
                        position.
  --PPG2KP-hybridize-with-restricted
                        If a cut matches the length or width of a
                        piece, and this length/width cannot be
                        obtained by any combination of smaller pieces,
                        then the cut already immediatelly extracts the
                        matching piece. Should preserve optimality but
                        it is experimental.
  --PPG2KP-faithful2furini2016
                        Tries to be the most faithful possible to the
                        description on the Furini2016 paper of the
                        PPG2KP (i.e., the model with the CutPosition
                        and RedundantCut reductions BUT NOT the
                        multistep pricing procedure); the flags
                        --no-cut-position, --no-redundant-cut,
                        --no-furini-symmbreak, can disable parts of
                        this reimplementation. Passing '--pricing
                        furini' will enable the pricing making the
                        model the Priced PPG2KP.
  --PPG2KP-no-redundant-cut
                        Disables Furini2016 Redundant-Cut reduction: a
                        bunch of flags is used to check if some trim
                        cut is necessary for optimality, or is
                        dominated by other trim cuts.
  --PPG2KP-no-cut-position
                        Disables Furini2016 Cut-Position reduction: if
                        a trivial heuristic shows that no combination
                        of six or more pieces fit a plate, then the
                        plate may be cut with restricted cuts without
                        loss to optimality.
  --PPG2KP-no-furini-symmbreak
                        Disables Furini2016 symmetry breaking: as trim
                        cuts are needed in faithful2furini2016 mode,
                        the symmetry-breaking is less restrictive than
                        'do not cut after midplate', it just removes
                        each cut y > midplate in which plate_size - y
                        is an existing cut; enabling this flag makes
                        the code just use all discretized positions as
                        cuts.
  --PPG2KP-pricing PPG2KP-PRICING
                        Define if a pricing procedure will be used,
                        and which one. The default value is 'none'.
                        The accepted values are: 'furini', 'becker',
                        and 'none'. Combining '--PPG2KP-pricing none'
                        with '--PPG2KP-faithful2furini2016' will give
                        the Complete PP-G2KP Model with Cut-Position
                        and Redundant-Cut reductions (pass
                        no-redundant-cut and no-cut-position to have
                        the pure Complete PP-G2KP Model). (default:
                        "none")
  --PPG2KP-ignore-d     Ignore the demand information during
                        discretization, used to measure impact. The
                        flag `hybridize-with-restricted` is affected
                        by it too (i.e., less hybridizations will
                        happen).
  --PPG2KP-use-c25      Add the tightening constraints 2.5. Increase
                        considerably the number of constraints hoping
                        to improve the relaxation a little.
  --PPG2KP-ignore-2th-dim
                        Ignore the dimension not being discretized
                        during discretization, used to measure impact.
                        The flag `hybridize-with-restricted` is
                        affected by it too (i.e., less hybridizations
                        will happen).
  --PPG2KP-heuristic-seed PPG2KP-HEURISTIC-SEED
                        Used only if the value of the 'MIP-start'
                        option implies the use of
                        `GuillotineModels.PPG2KP.Heuristic.fast_iterated_greedy`
                        (the seed is used to initialize a RNG, more
                        specifically a
                        `RandomNumbers.Xoroshiro128Plus`). (type:
                        Int64, default: 1)
  --PPG2KP-verbose      Only really relevant if you are studying the
                        code. Will increase the amount of data
                        outputted. The quiet flag overrides this.
  --PPG2KP-quiet        Avoid outputting any information or warning,
                        just output errors. Overrides the verbose
                        flag. .
  --PPG2KP-do-not-purge-unreachable
                        By default, the code removes variables and
                        constraints that are useless (because the only
                        way to reach them was removed by some other
                        mechanism, like pricing), this flag disables
                        this feature and leave the work to the
                        presolver of the MIP solver. Note that if the
                        'verbose' flag is enabled then the unreachable
                        variables and constraints will be searched and
                        printed, even if they are not removed from the
                        model after.
  --PPG2KP-MIP-start PPG2KP-MIP-START
                        Three values are allowed: "expected", "none",
                        and "guaranteed". If "expected" is passed,
                        then the model will not be MIP-started unless
                        some sort of pricing (that already generates a
                        complete feasible solution in the process) is
                        used. If it is Becker's pricing, then the
                        model is just MIP-started with the simple
                        iterated greedy 2-staged heuristic; if it is
                        Furini's pricing, then two MIP-starts occur:
                        the first one using the same heuristic in a
                        restricted model; the second one using the
                        solution of the restricted model in the
                        returned non-restricted model. If "none" is
                        passed, then the model is never MIP-started,
                        even temporarily (as the first MIP-start that
                        can occur in the middle of Furini's pricing).
                        If "guaranteed" is passed, it behaves as
                        "expected" if pricing is enabled, but if
                        pricing is disabled, it calls the iterated
                        greedy 2-staged heuristic just to MIP-start
                        the model before returning it. The MIP-start
                        is done by means of `MOI.set(model,
                        VariablePrimalStart(), var_index, value)`.
                        (default: "expected")
  --PPG2KP-allow-rotation
                        The formulation will allow the pieces to
                        rotate. Dummy rotated pieces are created and
                        feeded to the discretization and model
                        building. The only explicit difference in the
                        formulation is that the demand is shared
                        between the two rotations of a piece. Pieces
                        that meet any of the following conditions do
                        not have rotated counterparts: (i) the piece
                        is a square; (ii) after rotation the piece
                        does not fit into the 'original plate'/'large
                        object'; (iii) the instance has another piece
                        that is exactly the same as the rotated piece
                        (other attributes like profit need to match
                        too), in this case the demand of both
                        (non-rotated) pieces is summed and shared.
  --PPG2KP-mirror-plates
                        Will error if the allow-rotation flag is not
                        passed. If allow-rotation is enabled and
                        mirror-plates is disabled, then the cut/plate
                        enumeration is completely unaware of the fact
                        (pieces are duplicated before enumeration and
                        other changes are done at the formulation
                        level). If both flags are enabled, then the
                        enumeration is slightly changed: every plate
                        can only have length smaller than or equal to
                        the width (i.e., plates with length greater
                        than width are rotated automatocally). This
                        may cut by half the number of
                        plates/constraints and, consequently, also
                        reduce the number of cuts/variables (in other
                        words, the whole size of the model). The
                        Redundant-Cut improvement is disabled if this
                        flag has any effect (because we lose the
                        capability to distinguish between vertical and
                        horizontal cuts).
  --PPG2KP-aggressive-hybridization
                        Only can be used if
                        hybridization-with-restricted is enabled. The
                        default behavior of hybridization is to not
                        increase the model size significantly (a small
                        increase equivalent to the number of pieces
                        may occur). To guarantee this, however, we
                        need to abstain from hybridizations in a
                        specific circumnstance. The circumnstance is:
                        two or more pieces share the same length (or
                        width), and the associated discretization
                        point cannot be reached by a sum of smaller
                        pieces lengths (or widths). Without increasing
                        the number of variables (i.e., creating a
                        different hybridized cut for each piece), we
                        cannot hybridize for any of the pieces,
                        because hybridizing for a single arbitrary
                        piece would make the formulation lose the
                        guarantee of optimality. Consequently, the
                        default behaviour (i.e., without this flag)
                        loses some opportunities for hybridization to
                        keep optimality and to avoid considerable
                        model size increase. This flag allows the
                        formulation to seize these opportunities by
                        creating multiple hybridized cuts, it also
                        keeps the optimality but it can increase the
                        model size considerably.
```
