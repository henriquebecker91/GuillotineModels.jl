# = 2DPackGen =

# == Structures ==

# === Data/Problem Structures ===

struct SLOPP{D, S, P}
	"Length of the large object."
	L :: S
	"Width of the large object."
	W :: S
	"Lengths of the small objects."
	l :: Vector{S}
	"Widths of the small objects."
	w :: Vector{S}
	"Lower bound on piece demand (minimum amount needed for a valid solution)."
	dlb :: Vector{D}
	"Upper bound on piece demand (maximum amount allowed for a valid solution)."
	dub :: Vector{D}
	"Profits of the small objects; contribution to objective function if packed."
	p :: Vector{P}
end

struct MHLOPPW{D, S, P}
	"Lengths of the large objects."
	L :: Vector{S}
	"Widths of the large objects."
	W :: Vector{S}
	"Amount of copies available for each large object (cannot use more)."
	c :: Vector{S}
	"Lengths of the small objects."
	l :: Vector{S}
	"Widths of the small objects."
	w :: Vector{S}
	"Lower bound on piece demand (minimum amount needed for a valid solution)."
	lbd :: Vector{D}
	"Upper bound on piece demand (maximum amount allowed for a valid solution)."
	ubd :: Vector{D}
	"Profit of the pieces (contribution to objective function if packed)."
	p :: Vector{P}
end

struct ODPW{D, S, P}
	"Widths of the large objects."
	W :: Vector{S}
	"Amount of copies available for each large object (cannot use more)."
	copies :: Vector{S}
	"Cost associated with the use of CONTINUE HERE, CHECK WASCHER"
	cost :: Vector{S}
	"Lengths of the small objects."
	l :: Vector{S}
	"Widths of the small objects."
	w :: Vector{S}
	"Lower bound on piece demand (minimum amount needed for a valid solution)."
	lbd :: Vector{D}
	"Upper bound on piece demand (maximum amount allowed for a valid solution)."
	ubd :: Vector{D}
end

Number of different large objects (j)
LargeObject[j].Width	LargeObject[j].Available	LargeObject[j].Value
Number of different item types (i)
Item[i].Length	Item[i].Width	Item[i].Demand
# === Format/Instance Structures ===


function 


end

