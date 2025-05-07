using CSV
using DataFrames

function parse_load(df::DataFrame)
	# Parse the load data from the DataFrame and save it to CSV files
	# Load the required packages
	for (j,col) in enumerate([:P_load1, :P_load2, :P_load3, :P_load4])
		headers = []
		m = Vector{Vector{Float64}}(undef, length(df.lNodeID))
		for (i, row) in enumerate(eachrow(df))
	    		zone_id = row.zoneNum
	    		node_id = row.lNodeID
	    		header = string(col, "_", zone_id, "_", node_id)
	    		push!(headers, header)
	    		v = parse.(Float64, split(row[col][2:end-1], ","))
	    		m[i] = -v
		end
		CSV.write("/Users/sc87/code/OPF_LASCOPF_Staple/Horizontal_Proper/Example/load_data_stoc_$(col).csv", DataFrame(m, headers))
	end
end

load_df = CSV.File("/Users/sc87/code/OPF_LASCOPF_Staple/Horizontal_Proper/Example/load_data.csv") |> DataFrame
# Parse the load data
parse_load(load_df)
# Parse the load data from the DataFrame and save it to CSV files
# Load the required packages