struct TimeSection
	name :: String
	start :: Float64
end
export TimeSection

function Base.push!(
	stack :: Vector{TimeSection}, name :: String
)
	push!(stack, TimeSection(name, time()))
end

function Base.append!(
	stack :: Vector{TimeSection}, names :: AbstractVector{String}
)
	now = time()
	return append!(stack, (TimeSection(name, now) for name in names))
end

function close_and_print!(
	stack :: Vector{TimeSection}, name :: String,
	io = stdout
) :: Float64
	now = time()
	top = last(stack)
	@assert top.name == name
	println(io, "$(name) = $(now - top.start)")
	pop!(stack)
	return now
end
export close_and_print!

function close_and_print!(
	stack :: Vector{TimeSection},
	names :: Vector{String},
	io = stdout
) :: Float64
	now = time()
	to_pop_idxs = length(stack):-1:(length(stack) - length(names) + 1)
	for (idx_names, idx_stack) in enumerate(to_pop_idxs)
		if names[idx_names] != stack[idx_stack].name
			error(
				"$(names[idx_names]) != $(stack[idx_stack])." *
				"The sections must be closed in reverse order of opening (FIFO)."
			)
		end
	end
	for e in (stack[i] for i in to_pop_idxs)
		println(io, "$(e.name) = $(now - e.start)")
	end
	deleteat!(stack, reverse(to_pop_idxs)) # idxs must be sorted
	return now
end
export close_and_print!

