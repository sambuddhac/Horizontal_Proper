function load_input_files()
	shared_cand = DataFrame(XLSX.readtable("../../Input_Data/CandLine.csv", "Taul1", header=true)) #Dataframe of shared candidate lines
	int_cand =  DataFrame(XLSX.readtable("../../Input_Data/CandLineInt.csv", "Taul1", header=true)) #Dataframe of internal candidate lines
	shared_ex =  DataFrame(XLSX.readtable("../../Input_Data/SharedEline.csv", "Taul1", header=true)) #Dataframe of shared existing lines
	int_ex = DataFrame(XLSX.readtable("../../Input_Data/Tran.csv", "Taul1", header=true)) #Dataframe of internal existing lines
	gen =  DataFrame(XLSX.readtable("../../Input_Data/Gen.csv" , "Taul1", header=true)) #Dataframe of generators
	load =  DataFrame(XLSX.readtable("../../Input_Data/Load.csv", "Taul1", header=true)) #Dataframe of loads
	scen_prob = DataFrame(CSV.File("../../Input_Data/Scenario_Probability.csv", header=true)) #Dataframe of scenario probabilities
	zone_summary = DataFrame(CSV.File("../../Input_Data/Zone_Summary.csv", header=true)) #Dataframe of region-number of nodes

	input_dictionary = Dict()
	input_dictionary["shared_candidate_lines"] = shared_cand
	input_dictionary["internal_candidate_lines"] = int_cand
	input_dictionary["shared_existing_lines"] = shared_ex
	input_dictionary["internal_existing_lines"] = int_ex
	input_dictionary["generators"] = gen
	input_dictionary["load"] = load
	input_dictionary["scen_prob"] = scen_prob
	input_dictionary["zone_summary"] = zone_summary
end