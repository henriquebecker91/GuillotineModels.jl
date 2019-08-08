module Flow2Format

import LightGraphs, TikzGraphs, TikzPictures

push!(LOAD_PATH, "./")
using FlowDP

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
  num_nodes :: N,
  edges :: Vector{Edge{N, E}}
) :: TikzPictures.TikzPicture where {N, E}
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
    flow2digraph(num_nodes, edges), TikzGraphs.Layouts.SimpleNecklace(),
    options="scale=3, every node/.style={align=center}",
    edge_labels = d
  )
end

function flow2file(
  num_nodes :: N,
  edges :: Vector{Edge{N, E}},
  out :: TikzPictures.SaveType
) :: Nothing where {N, E}
  TikzPictures.save(out, flow2tikzpic(num_nodes, edges))
  nothing
end

function flow2pdf(
  num_nodes :: N,
  edges :: Vector{Edge{N, E}},
  fname :: String
) :: Nothing where {N, E}
  flow2file(num_nodes, edges, TikzPictures.PDF(fname))
  nothing
end

end # module

