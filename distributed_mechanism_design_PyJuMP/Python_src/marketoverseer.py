#Definition for Marketover class public Member Methods
import numpy as np
from Python_src.nettran import Nettran
from Python_src.log import log
from Python_src.profiler import Profiler
#define BIGM 1000000000000000000 // Defines the value of the Big M for transforming the bilinear terms to linear constraints

class Marketover(object):
	def __init__(self, zone_count, total_nodes, se_count, cand_count, vector_of_subnet, list_of_shared_node, list_of_global_shared_node, list_of_shared_zone, list_of_se_serial, list_of_se_from_rank, list_of_se_to_rank, list_of_se_react, list_of_se_cap, list_of_cand_serial, list_of_cand_from_rank, list_of_cand_to_rank, list_of_cand_react, list_of_cand_cap, milp_algo_choice, scenarios, compare_base): #constructor
		self.zone_number = zone_count #Initialize the total number of load zones
		self.sub_net_vector = vector_of_subnet #Create handle to the vector of subnetworks
		self.zone_list = list_of_shared_zone #List of zones between which exsiting and candidate shared lines exist
		self.node_list = list_of_shared_node #List of nodes at the ends of the existing and candidate shared lines
		self.shared_global_list = list_of_global_shared_node 
		self.node_number = total_nodes #total number of nodes, which form the ends of the shared lines
		self.cand_line_number = cand_count #Number of shared candidate lines
		self.shared_e_line_number = se_count #Number of shared existing Transmission lines
		self.se_serial = list_of_se_serial #Serial list of shared existing lines
		self.se_from_rank = list_of_se_from_rank #List of rank of from nodes of shared existing lines
		self.se_to_rank = list_of_se_to_rank #List of rank of to nodes of shared existing lines
		self.se_reactance = list_of_se_react #List of reactances of shared existing lines
		self.se_capacity = list_of_se_cap #List of line flow limits of shared existing lines
		self.cand_serial = list_of_cand_serial #Serial list of shared candidate lines
		self.cand_from_rank = list_of_cand_from_rank #List of rank of from nodes of shared candidate lines
		self.cand_to_rank = list_of_cand_to_rank #List of rank of to nodes of shared candidate lines
		self.cand_reactance = list_of_cand_react #List of reactances of shared candidate lines
		self.cand_capacity = list_of_cand_cap #List of line flow limits of shared candidate lines
		self.lp_solve_algo = milp_algo_choice #Simplex for 1 and IPM for 2
		self.count_of_scenarios = scenarios #Number of total random scenarios
		self.compare_basis = compare_base
		self.line_inter_dec_vec = np.zeros(int, float)
		self.ph_ang_dec_vec = np.zeros(int, float)
		"""
		int skip = 0;
		for (candserialiterator = cand_serial.begin(); candserialiterator != cand_serial.end(); ++candserialiterator) {
			if (skip > 0) {
				cout << "Candidate Serial is " << *candserialiterator << endl;
			}
			++skip;
		}
		for (candfromiterator = cand_from_rank.begin(); candfromiterator != cand_from_rank.end(); ++candfromiterator) {
			cout << "From rank is " << *candfromiterator << endl;
		}
		for (candtoiterator = cand_to_rank.begin(); candtoiterator != cand_to_rank.end(); ++candtoiterator) {
			cout << "To rank is " << *candtoiterator << endl;
		}
		"""
		#end constructor

	#def __del__(): #destructor
		#log.info("\nSimulation ended")
		#destructor ends

	def milp_market_overseer(self, lag_multi_xi, lag_multi_pi, total_cand_line_num, total_shared_node_num): #Function milp_market_overseer() implements the Mixed _integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
		#CREATION OF THE MIP SOLVER INSTANCE */
		dim_row = self.count_of_scenarios*(4 * self.cand_line_number + 2 * self.shared_e_line_number) #Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
		dim_col = self.count_of_scenarios*self.node_number+(self.count_of_scenarios+1)*self.cand_line_number #Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines

	#Function MILP() ends

	def lb_market_overseer(self, lag_multi_xi, lag_multi_pi, total_cand_line_num, total_shared_node_num): #Function lb_market_overseer() calculates the lower bound of the Mixed _integer Linear Programming Solver routine by calling GLPK routines for average heat rate objective
		dim_row = self.count_of_scenarios*(4 * self.cand_line_number + 2 * self.shared_e_line_number) #Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper line limits & lower and upper definition limits of candidate shared lines, and second term for lower and upper line limits for shared existing lines
		dim_col = self.count_of_scenarios*self.node_number+(self.count_of_scenarios+1)*self.cand_line_number #Total number of columns of the LP (number of Decision Variables) first term to account for voltage phase angles for inter-zonal lines' nodes, and second term for power flow values and binary integer decision variable values for shared candidate lines
	#Function MILP() ends

	def buffer_intermediate_decision(self, iter_count_out):
		if iter_count_out==0: 
			for scenario_track in range(self.count_of_scenarios): 
				for node_track in range(self.node_number):
					self.interim_cont_dec_var_prev.append(0)
	
			for line_track in range(self.cand_line_number):
				self.interim_int_dec_var_prev.append(0)
		else: 
			interim_cont_dec_var_prev = self.interim_cont_dec_var
			interim_int_dec_var_prev = self.interim_int_dec_var
	

	def get_global_upper(self, lag_multi_xi, lag_multi_pi, regional_upper, number_of_zones): #Function get_global_upper() returns the global upper bound for the investment coordination problem
		#Sum up the regional upper bounds, which are the tentative regional minima, at the end of every iteration
		revised_upper_bound = 0 #total upper bound initialized
		for i in range(number_of_zones):
			revised_upper_bound += regional_upper[i]
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
		for scen_pos in range(self.count_of_scenarios):
			rank_pres_checker = [] #Vector to check if a particular rank has been accounted for 
			length = len(self.angle_dec_index[scen_pos]) #length of the rank vector
			for i in range(length): #Put as many zeroes on the rank_pres_checker vector as is the length of the rank vector 
				rank_pres_checker.append(0)
			inter_angle_term_vec = [] #_intermediate vector for storing the costs of the angle terms from each sub-network
			tracker = 0 #tracker to track the position of the angle_dec_index_iterator
			for angle_dec_index_iterator in self.angle_dec_index[scen_pos]: #Iterate through rank vector
				inter_node_rank = angle_dec_index_iterator #Store the value of the rank of the present node iterate in the list 
				if rank_pres_checker[tracker] == 0: #If this node rank hasn't been already accounted for
					index_list = [pos1 for pos1 in range(len(self.angle_dec_index[scen_pos])) if self.angle_dec_index[scen_pos][pos1] == inter_node_rank] # Get all the indices of inter_node_rank in the angle_dec_index vector
					for index_list_iterator in index_list: # while all the different positions of this rank hasn't been accounted for
						rank_pres_checker[index_list_iterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
						inter_angle_term_vec.append(-(self.phase_angle_decision[scen_pos][index_list_iterator])*(lag_multi_xi[scen_pos*self.node_number+(angle_dec_index_iterator-1)])) #Calculate cost term
					smallest_element = min(inter_angle_term_vec)
					revised_upper_bound += smallest_element
					inter_angle_term_vec = []
					index_list = []
				tracker += 1 #Increment the tracker
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared discrete variable values from different regions/zones
		rank_pres_checker_int = [] #Vector to check if a particular rank has been accounted for 
		length_int = len(self.line_dec_index) #length of the rank vector
		for i in range(length_int):#Put as many zeroes on the rank_pres_checker vector as is the length of the rank vector 
			self.rank_pres_checker_int.append(0)
		inter_line_term_vec = [] #_intermediate vector for storing the costs of the line building decisions from each sub-network
		tracker_int = 0 #tracker to track the position of the line_dec_index_iterator
		for line_dec_index_iterator in self.line_dec_index: #Iterate through rank vector
			inter_line_rank = line_dec_index_iterator #Store the value of the rank of the present line iterate in the list 
			if rank_pres_checker_int[tracker_int] == 0: #If this line rank hasn't been already accounted for
				index_list = [pos1 for pos1 in range(len(self.line_dec_index)) if self.line_dec_index[pos1] == inter_line_rank] # Get all the indices of inter_node_rank in the angle_dec_index vector
				for index_list_iterator in index_list: # while all the different positions of this rank hasn't been accounted for
					rank_pres_checker_int[index_list_iterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					inter_line_term_vec.append(-(self.line_inter_decision[index_list_iterator])*(lag_multi_pi[(*line_dec_index_iterator-1)])) #Calculate cost term
				smallest_element = min(inter_line_term_vec)
				revised_upper_bound += smallest_element
				inter_line_term_vec = []
				index_list = []
			tracker_int += 1 #Increment the tracker
		return revised_upper_bound
	#Function get_global_upper ends

	def get_global_lower(regional_lower, number_of_zones): #Function get_global_lower() returns the global lower bound for the investment coordination problem
		revised_lower_bound = 0 #total lower bound initialized
		for i in range(number_of_zones):
			revised_lower_bound += regional_lower[i]
		return revised_lower_bound
	#Function get_global_lower ends

	def get_global_consensus(self): #Function get_global_consensus() returns the global consensus for the investment coordination problem
		glob_consensus = 0 #total consensus
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared continuous variable values from different regions/zones
		for scen_pos in range(self.count_of_scenarios):
			rank_pres_checker = [] #Vector to check if a particular rank has been accounted for 
			length = len(self.angle_dec_index[scen_pos]) #length of the rank vector
			for i in range(length): #Put as many zeroes on the rank_pres_checker vector as is the length of the rank vector 
				rank_pres_checker.append(0)
			tracker = 0 #tracker to track the position of the angle_dec_index_iterator
			for angle_dec_index_iterator in self.angle_dec_index[scen_pos]: #Iterate through rank vector
				inter_node_rank = angle_dec_index_iterator #Store the value of the rank of the present node iterate in the list 
				if rank_pres_checker[tracker] == 0: #If this node rank hasn't been already accounted for
					index_list = [pos1 for pos1 in range(len(self.angle_dec_index[scen_pos])) if self.angle_dec_index[scen_pos][pos1] == inter_node_rank] # Get all the indices of inter_node_rank in the angle_dec_index vector
					for index_list_iterator in index_list: # while all the different positions of this rank hasn't been accounted for
						rank_pres_checker[index_list_iterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
						glob_consensus += (self.phase_angle_decision[scen_pos][index_list_iterator]-self.interim_cont_dec_var[angle_dec_index_iterator-1]) ** 2 #Calculate cost term
					index_list = []
				tracker+=1 #Increment the tracker
		#Calculate the MO upper bound, which is the minimum value of the MO objective, calculated with a given mix of the intermediate, shared discrete variable values from different regions/zones
		rank_pres_checker_int = [] #Vector to check if a particular rank has been accounted for 
		length_int = len(self.line_dec_index) #length of the rank vector
		for i in range(length_int): #Put as many zeroes on the rank_pres_checker vector as is the length of the rank vector 
			rank_pres_checker_int.append(0)
		tracker_int = 0 #tracker to track the position of the line_dec_index_iterator
		for line_dec_index_iterator in self.line_dec_index: #Iterate through rank vector
			inter_line_rank = line_dec_index_iterator #Store the value of the rank of the present line iterate in the list 
			if rank_pres_checker_int[tracker_int] == 0: #If this line rank hasn't been already accounted for
				index_list = [pos1 for pos1 in range(len(self.line_dec_index)) if self.line_dec_index[pos1] == inter_line_rank] # Get all the indices of inter_node_rank in the angle_dec_index vector
				for index_list_iterator in index_list: # while all the different positions of this rank hasn't been accounted for
					rank_pres_checker_int[index_list_iterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					glob_consensus += (self.line_inter_decision[index_list_iterator]-self.interim_int_dec_var[line_dec_index_iterator-1]) ** 2 #Calculate cost term
				index_list = []
			tracker_int += 1 #Increment the tracker
		return glob_consensus**0.5
	
	#Function get_global_consensus ends

	def finDecLineConstr(self): #Final decisions on the construction stauses of candidate lines taking into account the decisions of different zones
		rank_pres_checker_int = [] #Vector to check if a particular rank has been accounted for 
		length_int = len(self.line_dec_index) #length of the rank vector
		for i in range(length_int): #Put as many zeroes on the rank_pres_checker vector as is the length of the rank vector 
			rank_pres_checker_int.append(0)
		tracker_int = 0 #tracker to track the position of the line_dec_index_iterator
		for line_dec_index_iterator in self.line_dec_index: #Iterate through rank vector
			inter_line_rank = line_dec_index_iterator #Store the value of the rank of the present line iterate in the list 
			if rank_pres_checker_int[tracker_int] == 0: #If this line rank hasn't been already accounted for
				first = 0 #Initialize a flag to indicate the first occurence of a rank
				first_in_series = 0 #flag to indicate if the first verdict for this rank is 1 or 0
				index_list = [pos1 for pos1 in range(len(self.line_dec_index)) if self.line_dec_index[pos1] == inter_line_rank] # Get all the indices of inter_node_rank in the angle_dec_index vector
				for index_list_iterator in index_list: # while all the different positions of this rank hasn't been accounted for
					rank_pres_checker_int[index_list_iterator] = 1 #Mark the element in the corresponding position of checker vector 1, to indicate this rank has been accounted for 
					if self.line_inter_decision[index_list_iterator]==1 and first == 0: #Check if the decision of the different subnets is 1
						self.constructed_ranks.append(inter_line_rank)
						first_in_series = 1 #Set the flag to indicate that the first verdict for this rank is 1
					elif self.line_inter_decision[index_list_iterator]==0 and first != 0 and first_in_series == 1:
						self.constructed_ranks.pop()
					first += 1 #Increment the flag to indicate the recurrence of the rank
				index_list = []
			tracker_int += 1 #Increment the tracker
		"""
		for construc_iterator in self.constructed_ranks:	
			log.info(" First message ")
			log.info("Constructed shared candidate lines are {}".format(construc_iterator))
			log.info(" Second message ")
		log.info("end of finDecLineConstr ")
		"""

	def scan_built_lines_list(self, rank_glob): #Scans the list of built candidate lines for the input argument of global rank to check if this line is built
		pos = self.constructed_ranks.index(rank_glob) #find the first position of this rank in the vector
		if pos != self.constructed_ranks[len(self.constructed_ranks)-1]: #if the rank is found in the list
    			return 1 #return 1
		else:
			return 0 #Otherwise, return 0

	def populate_line_dec(self, interim_dec_bin, zonal_index, global_rank):#Method to pass the intermediate message for the zonal line building decision to the MO
		self.line_inter_decision.append(interim_dec_bin)
		self.line_dec_index.append(global_rank)
		self.zonal_ind_vector_int.append(zonal_index)
		self.line_inter_dec_vec[zonal_index*self.cand_line_number+global_rank]=interim_dec_bin

	def populate_angle_dec(self, interim_dec_angle, zonal_index, scenario, global_rank): #Method to pass the intermediate message for the zonal node angle decision to the MO
		(self.phase_angle_decision[scenario]).append(interim_dec_angle)
		(self.angle_dec_index[scenario]).append(global_rank)
		(self.zonal_ind_vector_cont[scenario]).append(zonal_index)
		self.ph_ang_dec_vec[scenario*self.zone_number*self.node_number+zonal_index*self.node_number+global_rank]=interim_dec_angle

	def reward_penalty_cont(self, lagrange_mult, matrix_index, iteration): #Function get_global_upper() returns the global upper bound for the investment coordination problem
		#cout << "\nITERATION COUNT : " << iteration << endl;
		for scen_pos in range(self.count_of_scenarios):
			tracker = 0
			for angle_dec_index_iterator in self.angle_dec_index[scen_pos]:
				if matrix_index==(scen_pos*self.zone_number*self.node_number+(self.zonal_ind_vector_cont[scen_pos][tracker]-1)*self.node_number+(angle_dec_index_iterator)):		
					if self.compare_basis == 1:
						lagrange_mult += (self.phase_angle_decision[scen_pos][tracker]-self.interim_cont_dec_var[scen_pos*self.node_number+angle_dec_index_iterator-1])
						#cout << "From MO Continuous Lagrange Multiplier update: " << "Angle from Zones. Tracker: " << tracker << " Scenario: " << scen_pos << " Zonal Angle value " << ((phase_angle_decision[scen_pos]).at(tracker)) << " MO Angle index: " << scen_pos*node_number+(*angle_dec_index_iterator-1) << " MO Angle value: " << (interim_cont_dec_var[scen_pos*node_number+(*angle_dec_index_iterator-1)]) << " Updated Lagrange Multiplier: " << lagrange_mult << endl;
					else:
						lagrange_mult += (self.phase_angle_decision[scen_pos][tracker]-self.interim_cont_dec_var_prev[scen_pos*self.node_number+angle_dec_index_iterator-1])
						#log.info("From MO Continuous Lagrange Multiplier update: Angle from Zones. Tracker: {} Scenario: {} Zonal Angle value {} MO Angle index: {} MO Angle value: {} Updated Lagrange Multiplier: {}".format(tracker, scen_pos, self.phase_angle_decision[scen_pos][tracker], scen_pos*self.node_number+(angle_dec_index_iterator-1), interim_cont_dec_var[scen_pos*node_number+(angle_dec_index_iterator-1)], lagrange_mult))
				tracker += 1
			#log.info("Tracker {}".format(tracker))
		return lagrange_mult
	#Function get_global_upper ends

	def reward_penalty_integer(self, lagrange_mult, matrix_index, iteration): #Function get_global_upper() returns the global upper bound for the investment coordination problem
		#log.info("\nITERATION COUNT : {}".format(iteration))
		tracker = 0
		for line_dec_index_iterator in self.line_dec_index:
			if matrix_index==(self.zonal_ind_vector_int[tracker]-1)*self.cand_line_number+line_dec_index_iterator:
				if self.compare_basis == 1:	
					lagrange_mult += (self.line_inter_decision[tracker]-self.interim_int_dec_var[line_dec_index_iterator-1])
					#cout << "From MO _integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << line_inter_decision[tracker] << " MO decision index: " << (*line_dec_index_iterator-1) << " MO Decision value: " << interim_int_dec_var[(*line_dec_index_iterator-1)] << " Updated Lagrange Multiplier: " << lagrange_mult << endl
				else:
					lagrange_mult += (self.line_inter_decision[tracker]-self.interim_int_dec_var_prev[line_dec_index_iterator-1])
					#log.info("From MO _integer Lagrange Multiplier update: " << "Line Decision from Zones. Tracker: " << tracker << " Zonal Decision value " << line_inter_decision[tracker] << " MO decision index: " << (*line_dec_index_iterator-1) << " MO Decision value: " << interim_int_dec_var[(*line_dec_index_iterator-1)] << " Updated Lagrange Multiplier: " << lagrange_mult << endl
			tracker+=1
		#log.info("Tracker {}".format(tracker))
		return lagrange_mult	
	#Function get_global_upper ends

	def clear_vectors(self): #Clears the different interim vectors for making them ready for the next iteration
		for scen_pos in range(self.count_of_scenarios):
			#Clear vectors for interim continuous decision variables
			self.phase_angle_decision[scen_pos] = []
			self.angle_dec_index[scen_pos] = []
			self.zonal_ind_vector_cont[scen_pos] = []
		self.interim_cont_dec_var = []
		self.interim_cont_dec_var_prev = []
		#Clear vectors for interim integer decision variables
		self.line_inter_decision = []
		self.line_dec_index = []
		self.zonal_ind_vector_int = []
		self.interim_int_dec_var = []
		self.interim_int_dec_var_prev = []

	def clear_delayed_vectors(self): #Clears the different interim vectors only buffer vectors
		self.interim_cont_dec_var_prev = []
		self.interim_int_dec_var_prev = []
