module Format

# CONSIDER MIGRATING TO GraphRecipes
# (https://github.com/JuliaPlots/GraphRecipes.jl/blob/master/README.md)
# pros: seems more mature, use the Plots library
# cons: seems even less documented

import LightGraphs, TikzGraphs, TikzPictures

using ..Enumeration

export flow2digraph, flow2tikzpic, flow2file, flow2pdf

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

end # module

