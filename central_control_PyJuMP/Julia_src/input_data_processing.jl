function input_data_processing(indict::Dict)
	shared_cand = indict["shared_candidate_lines"]
	int_cand = indict["internal_candidate_lines"]
	shared_ex = indict["shared_existing_lines"]
	int_ex = indict["internal_existing_lines"]
	gen = indict["generators"]
	load = indict["load"]
	scen_prob = indict["scen_prob"]
	
	l(i,s) = load[load.zoneNum .== i, [1,2,s+2]] # load within zone i
	g(i) = gen[gen.zoneNum .== i, :]   # generators within zone i
	shared_c(i) = vcat(shared_cand[shared_cand.nodeZone1 .== i,:] , shared_cand[shared_cand.nodeZone2 .== i, :]) #shared candidate lines within zone i
	int_c(i) = int_cand[int_cand.zoneNum .== i, :]   # number of internal candidate lines within zone i
	shared_e(i) = vcat(shared_ex[shared_ex.nodeZone1 .== i,:] , shared_ex[shared_ex.nodeZone2 .== i, :]) #shared existing lines within zone i
	int_e(i) =int_ex[int_ex.zoneNum .== i, :]       # internal existing lines within zone i
	MC(i) = (g(i).C2 .* (g(i).PgMax .^ 2) .+ g(i).C1 .* g(i).PgMax .- g(i).C2 .*(g(i).PgMin .^ 2) .- g(i).C1.* g(i).PgMin) ./ (g(i).PgMax .- g(i).PgMin) #Marginal cost of generators within zone i
	bin_c(i) = (shared_cand.nodeZone1 .== i) + (shared_cand.nodeZone2 .== i) # A binary vector through which we can check if the shared candidate lines belong to zone i
	bin_e(i) = (shared_ex.nodeZone1 .== i) + (shared_ex.nodeZone2 .== i) # A binary vector through which we can check if the shared existing lines belong to zone i
	scen_weight(s) = scen_prob.scen_weight[scen_prob.scenario .== s,:]
end