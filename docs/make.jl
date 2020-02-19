using Documenter, GuillotineModels

makedocs(
  modules = [GuillotineModels]
  , format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"
  )
  , checkdocs = :exports
  , sitename = "GuillotineModels.jl"
  , pages = Any[
		"GuillotineModels" => "index.md",
		"GuillotineModels.Utilities" => [
			"Utilities" => "Utilities.md",
			"Utilities.Args" => "U-Args.md"
		],
		"GuillotineModels.InstanceReader" => "InstanceReader.md",
		"GuillotineModels.CommandLine" => [
			"CommandLine" => "CommandLine.md",
			"CommandLine.SolversArgs" => "CL-SolversArgs.md"
		],
		"GuillotineModels.PPG2KP" => [
			"PPG2KP" => "PPG2KP.md",
			"PPG2KP.Heuristic" => "PPG2KP-Heuristic.md",
			"PPG2KP.Enumeration" => "PPG2KP-Enumeration.md"
		],
		"GuillotineModels.Flow" => [
			"Flow" => "Flow.md",
			"Flow.Enumeration" => "Flow-Enumeration.md",
			"Flow.Format" => "Flow-Format.md"
		]
	]
  , doctest = true
)

deploydocs(
  repo = "github.com/henriquebecker91/GuillotineModels.jl.git",
)
