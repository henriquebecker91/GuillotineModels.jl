const _V_OPEN = "["
const _H_OPEN = "{"
_is_open(s) = s == _V_OPEN || s == _H_OPEN

# We have a non-capturing and a capturing version of each regex.
# The non-capturing are used for checking validity, the capturing are
# used for parsing. Ideally we should only use the capturing version
# and cache the result but this is not performance-sensitive code.
# So our trade-off here is just not doing the capturing (slower) version
# two times for each token.
const _PIECE_REGEX_NC = r"^[1-9]\d*p[1-9]\d*x[1-9]\d*$"
const _PIECE_REGEX_C = r"^([1-9]\d*)p([1-9]\d*)x([1-9]\d*)$"
const _PLATE_REGEX_NC = Regex(
	"^P[1-9]\\d*x[1-9]\\d*[$(_V_OPEN)$(_H_OPEN)]?\$"
)
const _PLATE_REGEX_C = Regex(
	"^P([1-9]\\d*)x([1-9]\\d*)([$(_V_OPEN)$(_H_OPEN)]?)\$"
)

const _V_CLOSE = "]"
const _H_CLOSE = "}"
_is_close(s) = s == _V_CLOSE || s == _H_CLOSE

_is_valid_piece(s) = occursin(_PIECE_REGEX_NC, s)
_is_valid_plate(s) = occursin(_PLATE_REGEX_NC, s)
# The valid tokens do not include _is_open because these characters are
# part of a valid plate token.
_is_valid_token(s) = _is_close(s) || _is_valid_piece(s) || _is_valid_plate(s)

# Only guarantee correctness and not throwing an exception if the token
# returns true for `_is_valid_token`
_is_piece(token) = isdigit(first(token))
_is_plate(token) = first(token) === 'P'

function _parse_piece(token, ::Type{D}, ::Type{S}) where {D, S}
	m = match(_PIECE_REGEX_C, token)
	i, l, w =
		parse(D, m.captures[1]), parse(S, m.captures[2]), parse(S, m.captures[3])

	return CutPattern(l, w, i)
end

function _are_matching_brackets(open, close)
	open == "[" && close == "]" && return true
	open == "{" && close == "}" && return true
	return false
end

function _parse_plate!(
	token, pattern_stack :: Vector{Vector{CutPattern{D, S}}}, close_stack
) where {D, S}
	m = match(_PLATE_REGEX_C, token)
	L, W = parse(S, m.captures[1]), parse(S, m.captures[2])
	open = m.captures[3]
	if isempty(open) # i.e., the plate is waste
		subpatterns = CutPattern{D, S}[]
		is_vertical = true
	else # i.e., the plate has subpatterns
		subpatterns = pop!(pattern_stack)
		@assert _is_open(open)
		close = pop!(close_stack)
		if !_are_matching_brackets(open, close)
			error(
				"A subpattern started by a $(open) bracket ends with a" *
				" $(close) bracket.")
		end
		is_vertical = open == _V_OPEN
	end

	return CutPattern(L, W, is_vertical, subpatterns)
end

function parse_pattern(
	str :: String, ::Type{D}, ::Type{S}
) :: Vector{CutPattern{D, S}} where {D, S}
	tokens = split(str)

	# Note: pattern_stack starts with an empty level, while close_stack starts
	# completely empty.
	pattern_stack = Vector{CutPattern{D, S}}[CutPattern{D, S}[]]
	close_stack = String[]

	for token in reverse(tokens)
		if !_is_valid_token(token)
			error("Token \"$token\" was not recognized as a valid token.")
		end

		if _is_close(token)
			push!(pattern_stack, CutPattern{D, S}[])
			push!(close_stack, token)
			continue
		end

		if _is_piece(token)
			push!(last(pattern_stack), _parse_piece(token, D, S))
			continue
		end

		@assert _is_plate(token)

		# A plate may have subpatterns or just be waste. If it is waste,
		# then pattern_stack and close_stack are left untouched. If it has
		# subpatterns then it check the correctness by popping close_stack
		# and it gets its subpatterns by popping pattern_stack.
		plate = _parse_plate!(token, pattern_stack, close_stack)
		push!(last(pattern_stack), plate)
	end

	if !isempty(close_stack)
		error("The pattern finished without closing all subpatterns.")
	end

	return only(pattern_stack)
end

