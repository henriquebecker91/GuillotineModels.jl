using Documenter, GuillotineModels

makedocs(
  modules = [GuillotineModels]
  , format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"
  )
  , checkdocs = :exports
  , sitename = "GuillotineModels.jl"
  , pages = Any[
		"GuillotineModels.md",
		"Utilities" => [
			"Utilities.md",
			"U-Args.md"
		],
		"InstanceReader" => "InstanceReader.md",
		"CommandLine" => [
			"CommandLine.md",
			"CommandLine.SolverArgs" => "CL-SolverArgs.md"
		],
		"PPG2KP" => [
			"PPG2KP.md",
			"PPG2KP.Heuristic" => "PPG2KP-Heuristic.md",
			"PPG2KP.Enumeration" => "PPG2KP-Enumeration.md"
		],
		"Flow" => [
			"Flow.md",
			"Flow.Enumeration" => "Flow-Enumeration.md",
			"Flow.Format" => "Flow-Format.md"
		]
	]
  , doctest = true
)

deploydocs(
  repo = "github.com/henriquebecker91/GuillotineModels.jl.git",
)
