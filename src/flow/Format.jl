"""
Flow.Format consists uniquely of optional methods for making it easier to
visualize the structure of the Flow model. It requires the user to install and
import packages not in the GuillotineModels dependencies in order to have
access to its methods. The `flow2digraph` method needs LightGraphs, and the
other three (`flow2file`, `flow2tikzpic`, `flow2pdf`) need LightGraphs,
TikzGraphs, and TikzGraphs.

The `Node` and `Edge` types are defined by Flow.Enumeration and the Vectors
of all nodes and all edges of the model is returned by
`GuillotineModels.build_model(::Val{:Flow}, ...)`.
"""
module Format

# CONSIDER MIGRATING TO GraphRecipes
# (https://github.com/JuliaPlots/GraphRecipes.jl/blob/master/README.md)
# pros: seems more mature, use the Plots library
# cons: seems even less documented

using ..Enumeration

import Requires.@require

function __init__()
	@require LightGraphs="093fc24a-ae57-5d10-9952-331d41423f4d" begin
		export flow2digraph
		function flow2digraph(
			num_nodes :: N,
			edges :: Vector{Edge{N, E}}
		) :: LightGraphs.SimpleDiGraph{E} where {N, E}
			g = LightGraphs.SimpleDiGraph{E}(num_nodes)

			for e in edges
				println(e.indx, ": ", e.head, " -> ", e.tail)
				LightGraphs.add_edge!(g, e.head, e.tail)
			end

			g
		end

		@require TikzPictures="37f6aa50-8035-52d0-81c2-5a1d08754b2d" begin
		@require TikzGraphs="b4f28e30-c73f-5eaf-a395-8a9db949a742" begin
			export flow2tikzpic, flow2file, flow2pdf
			function flow2tikzpic(
				nodes :: Vector{Node{S, N}},
				edges :: Vector{Edge{N, E}}
			) :: TikzPictures.TikzPicture where {S, N, E}
				@assert issorted(nodes, by = n -> n.idx)
				d = Dict{Tuple{E, E}, String}()
				for e in edges
					arc = convert.(E, (e.head, e.tail))
					if !haskey(d, arc)
						d[arc] = "$(e.indx)+$(e.back)"
					else
						d[arc] *= ", \\\\ $(e.indx)+$(e.back)"
					end
				end

				TikzGraphs.plot(
					flow2digraph(length(nodes), edges),
					labels = map(n -> "$(n.idx):$(n.par):$(n.per):$(n.ori)", nodes),
					TikzGraphs.Layouts.SimpleNecklace(),
					options="scale=3, every node/.style={align=center}",
					edge_labels = d
				)
			end

			function flow2file(
				nodes :: Vector{Node{S, N}},
				edges :: Vector{Edge{N, E}},
				out :: TikzPictures.SaveType
			) :: Nothing where {S, N, E}
				TikzPictures.save(out, flow2tikzpic(nodes, edges))
				nothing
			end

			function flow2pdf(
				nodes :: Vector{Node{S, N}},
				edges :: Vector{Edge{N, E}},
				fname :: String
			) :: Nothing where {S, N, E}
				flow2file(nodes, edges, TikzPictures.PDF(fname))
				nothing
			end
		end # TikzPictures
		end # TikzGraphs
	end # LightGraphs
end # __init__
end # module

