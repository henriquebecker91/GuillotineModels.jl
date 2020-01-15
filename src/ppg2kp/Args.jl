module Args

function get_arg_parse_settings()
	s = ArgParseSettings()
	@add_arg_table s begin
	end
	s
end

end # module
