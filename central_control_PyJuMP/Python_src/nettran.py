#Definition for Nettran class
import julia
import os
import subprocess
import pandas as pd
import numpy as np
import json
import sys
import traceback
from Python_src.log import log
from Python_src.profiler import profiler
import gurobipy as gp
from gurobipy import GRB
from Python_src.powergenerator import power_generator
from Python_src.transl import transmission_line
from Python_src.load import Load
from Python_src.node import Node
from Python_src.sharedLine import se_line
from Python_src.candidateLine import cand_line
from Python_src.intcandidateLine import int_cand_line
#define AVERAGE_HEAT 1 // Defines the Average Heat generator cost function mode
#define PIECEWISE_LINEAR 2 // Defines the Piecewise Linear generator cost function mode
#define POLYNOMIAL 3 // Defines the Convex Polynomial generator cost function mode
#define BIGM 1000000000000000000 // Defines the value of the Big M for transforming the bilinear terms to linear constraints

profiler = profiler()

if sys.platform in ["darwin", "linux"]:
	log.info("Using Julia executable in {}".format(str(subprocess.check_output(["which", "julia"]), 'utf-8').strip('\n')))
elif sys.platform in ["win32", "win64", "cygwin"]:
	log.info("Using Julia executable in {}".format(str(subprocess.check_output(["which", "julia"]), 'utf-8').strip('\n')))

log.info("Loading Julia...")
profiler.start()
julSol = julia.Julia()
julSol.using("Pkg")
julSol.eval('Pkg.activate(".")')
julSol.include(os.path.join("JuMP_src", "HorMILPDistMech.jl")) # definition of Gensolver class for base case scenario first interval
log.info(("Julia took {:..2f} seconds to start and include Horizontal Investment Coordination mechanism design models.".format(profiler.get_interval())))

class Nettran(object):
	def __init__(self, json_list, zone_count, obj_choice): # constructor
		self.sim_mode = obj_choice #Specify the type of the curve
		self.zonal_count = zone_count
		self.node_num_vector = []
		self.node_num_vector.append(0) #Initialize the node number vector to indicate there no nodes in the fictitous 0-th zone
		for zonal_index in json_list: #Iterate through the zones 
			self.net_file = zonal_index['Network File'] #String for storing the name of the network file
			self.gen_file = zonal_index['Generator File'] #String for storing the name of the generator file
			self.tran_file = zonal_index['Transmission Lines File'] #String for storing the name of the transmission line file
			self.load_file = zonal_index['Load File'] #String for storing the name of the load file
			self.int_cand_line_file = zonal_index['Intra Candidate Lines File'] #String for storing the name of the candidate lines file
			#/* Nodes */
			matrix_net_file = json.load(open(os.path.join("data", self.net_file)))
			self.node_number = matrix_net_file['node_number']
			self.gen_number = matrix_net_file['gen_number'] 
			self.load_number = matrix_net_file['load_number'] 
			self.tran_number = matrix_net_file['tran_number']
			self.count_of_scenarios = matrix_net_file['load_stoch_scenarios'] #get the number of stochastic scenarios of load#get the dimensions of the Network
			for l in range(self.node_number):
				#log.info("Creating the {} -th Node".format(l + 1)) 
		
				node_instance = Node( self.univ_node_num + l + 1, l + 1, zonal_index ) #creates node_instance object with ID l + 1

				self.node_object.append( node_instance ) #pushes the node_instance object into the vector
			#end initialization for Nodes
##/* Generators */
			#/* Instantiate Generators */
        	matrix_gen = json.load(open(os.path.join("data", self.gen_file))) #ifstream constructor opens the file of Generators
			j = 0 #counter for generators
			for matrix_gen_file in matrix_gen:
				g_node_id = matrix_gen_file['gen_node_id'] #node object ID to which the particular generator object is connected
			#log.info("\nConnection Node defined.\n")
			#Parameters for Generator
			#Quadratic Coefficient: 
				c2 = matrix_gen_file['quad_cost_coeff'] * (100**2)
			#Linear coefficient: 
				c1 = matrix_gen_file['lin_cost_coeff'] * 100
			#Constant term: 
				c0 = matrix_gen_file['no_load_cost']
			#Maximum Limit: 
				pg_max= matrix_gen_file['pg_max'] / 100
			#Minimum Limit: 
				pg_min = matrix_gen_file['pg_min'] / 100
			#/* Secant Approximation of the Quadratic Cost curve */
			#Tangent Ratio of the secant approximation of the intercepted cost curve
				tan_theta = (c2*(pg_max**2)+c1*pg_max-c2*(pg_min**2)-c1*pg_min)/(pg_max-pg_min)
				min_cost = c2*(pg_min**2)+c1*pg_min+c0-tan_theta*pg_min
			#Intercept value or cost at minimum power level
			#check the bounds and validity of the parameter values
				gen_instance = self.power_generator( univ_gen_num + j+1, node_object[ univ_node_num + g_node_id - 1 ],  tan_theta, min_cost, pg_max, pg_min )
				self.gen_object.append(gen_instance) #push the generator object into the array
				j +=1 #increment counter
			matrix_gen.close() #Close the generator file
			self.gen_df = pd.dataframe([g.__dict__ for g in self.gen_object])
##/* Transmission Lines */
        		if self.tran_number > 0:
				matrix_tran_file = json.load(open(os.path.join("data", self.tran_file))) #ifstream constructor opens the file of Transmission lines
		#/* Instantiate Transmission Lines */
				k=0 #counter for transmission lines
				for matrix_tran in matrix_tran_file:
			#node IDs of the node objects to which this transmission line is connected.
					t_node_id1 = matrix_tran['from_node'] #From end
					t_node_id2 = matrix_tran['to_node'] #To end
			#reactance
					react = matrix_tran['reactance']
			#values of maximum allowable power flow on line in the forward and reverse direction:
					pt_max = matrix_tran['line_limit']/100
			#creates trans_line_instance object with ID k + 1
					trans_line_instance = self.transmission_line( univ_tran_num + k + 1, self.node_object[ univ_node_num + t_node_id1 - 1 ], self.node_object[ univ_node_num + t_node_id2 - 1 ], pt_max, react ) 
					self.transl_object.append( trans_line_instance ) #pushes the trans_line_instance object into the vector
					k +=1 #increment the counter
		#end initialization for Transmission Lines 
				matrix_tran_file.close() #Close the transmission line file 
				self.tran_df = pd.dataframe([t.__dict__ for t in self.transl_object])
##/* Internal Candidate Transmission Lines */
			if self.internal_c_lines > 0:
        			matrix_int_ce_tran_file = json.load(open(os.path.join("data", self.int_cand_line_file))) #ifstream constructor opens the file of internal candidate Transmission lines
        			k = 0 #Counter for internal candidate lines
			#/* Instantiate Internal Candidate Transmission Lines */
				for matrix_int_ce_tran in matrix_int_ce_tran_file:
				#node object IDs to which the particular transmission line object is connected
				#node IDs of the node objects to which this transmission line is connected.
					t_node_id1 = matrix_int_ce_tran['from_node'] #From end node 
					t_node_id2 = matrix_int_ce_tran['to_node'] #To end node
				#Parameters for Transmission Line
				#reactance:
					react = matrix_int_ce_tran['reactance']
				#values of maximum allowable power flow on line in the forward and reverse direction:
				#Forward direction:
					pt_max = matrix_int_ce_tran['line_limit']/100
					life_time = matrix_int_ce_tran['life_time'] #life time of the candidate line
					interest_rate = matrix_int_ce_tran['interest_rate'] #interest rate of the investment 
					cost_per_cap = matrix_int_ce_tran['cost_per_cap']*pt_max #capital cost for the construction 
					pres_absence = matrix_int_ce_tran['pres_absence'] #status of the construction
			
				#creates int_cand_line_instance object with ID k + 1
					int_cand_line_instance = self.int_cand_line( univ_int_cand_num + k + 1, self.node_object[ univ_node_num + t_node_id1 - 1 ], self.node_object[ univ_node_num + t_node_id2 - 1 ], pt_max, react, interest_rate, life_time, cost_per_cap, pres_absence ) #Create the internal candidate transmission line object with node 1 
					self.int_cand_line_object.append( int_cand_line_instance ) #pushes the trans_line_instance object into the vector
            				k +=1 #increment the counter
			#end initialization for candidate Transmission Lines
				matrix_int_ce_tran_file.close() #Close the candidate lines file
				self.int_cand_tran_df = pd.dataframe([t.__dict__ for t in elf.int_cand_line_object])
##/* Loads */
        		matrix_load_file = json.load(open(os.path.join("data", self.load_file))) #ifstream constructor opens the file of Loads
			count_of_scenarios = load_fields-1
			for l in range(self.node_number+self.node_number):
		    		(self.node_object[l]).int_load(self.count_of_scenarios) #Initialize the default loads on all nodes to zero
		#end initialization for Nodes		
		#/* Instantiate Loads */
			for j in load_number:
				for matrix_load in matrix_load_file:
			#node object ID to which the particular load object is connected
			#node ID of the node object to which this load object is connected.
					l_node_id = matrix_load['load_node_id']
					for f in range(self.count_of_scenarios):
				#value of allowable power consumption capability of load in pu with a negative sign to indicate consumption:
				#Power Consumption:
						p_load[f] = matrix_load['load_mw']/100

				load_instance = self.Load( univ_load_num + j + 1, self.node_object[ univ_node_num + l_node_id - 1 ], load_fields-1, p_load ) #creates load_instance object object with ID number j + 1

				self.load_object.append( load_instance ) #pushes the load_instance object into the vector
		#end initialization for Loads
			matrix_load_file.close() #Closes the load file
			for f in range(self.count_of_scenarios):
				self.probability.append(1/(self.count_of_scenarios))
########################################################################################################################	
		univ_node_num = univ_node_num + node_number #Increment the universal node number
		node_num_vector.append(node_number)
		univ_gen_num = univ_gen_num + gen_number #Increment the universal generator number
		univ_tran_num = univ_tran_num + tran_number #Increment the universal transmission line number
		univ_load_num = univ_load_num + load_number #Increment the universal load number
		univ_int_cand_num = univ_int_cand_num + internal_c_lines #Increment the universal intra zonal candidate line number
######################################################################################################################## 
		for zonal_index in json_list: #Iterate through the zones 
			self.net_file = zonal_index['Network File'] 
			self.shared_line_file = json_index['Shared Lines File'] 
			self.tran_file = zonal_index['Transmission Lines File'] 
			self.cand_line_file  = json_index['Candidate Lines File'] 
## Shared Existing Transmission lines
			matrix_se_tran_file = json.load(open(os.path.join("data", self.shared_line_file))) #ifstream constructor opens the file of Transmission lines
			k = 0
		#/* Instantiate Shared Existing Transmission Lines */
			for matrix_se_tran in matrix_se_tran_file:
			#log.info("Tran File Test Message 1 from line {} before creation1".format(k))
			#node IDs of the node objects to which this transmission line is connected.
				ser_num = matrix_se_tran['global_serial'] #global serial number of the shared existing transmission line
				t_node_id1 = matrix_se_tran['from_node'] #From end node 
				node_zone1 = matrix_se_tran['from_zone'] #From end zone number
				t_node_id2 = matrix_se_tran['to_node'] #To end node
				node_zone2 = matrix_se_tran['to_zone'] #To end zone number
			#reactance:
				react = matrix_se_tran['reactance']
			#values of maximum allowable power flow on line in the forward and reverse direction:
				pt_max = matrix_se_tran['line_limit']/100
				if ser_num in se_line_ser_list:
					new_node_base = 0
					for aggr_count in node_zone1:
						new_node_base = new_node_base + node_num_vector[aggr_count]
					t_node_id1 = t_node_id1 + new_node_base
					new_node_base = 0
					for aggr_count in node_zone2:
						new_node_base = new_node_base + node_num_vector[aggr_count]
					t_node_id2 = t_node_id2 + new_node_base
					se_line_instance = self.se_line( k + 1, ser_num, self.node_object[ t_node_id1 - 1 ], self.node_object[ t_node_id2 - 1 ], pt_max, react )
					se_line_objectt.append( se_line_instance )
					se_line_ser_list.append(ser_num)
					univse_line_num = univse_line_num + 1
					k +=1
			matrix_se_tran_file.close() #Close the shared existing file
			self.shared_existing_df = pd.dataframe([se_df.__dict__ for se_df in self.se_line_objectt])
##/* Shared Candidate Transmission Lines */
        		matrix_ce_tran_file = json.load(open(os.path.join("data", self.cand_line_file ))) #ifstream constructor opens the file of candidate Transmission lines
		#/* Instantiate Shared Candidate Transmission Lines */
			k=0
			for matrix_ce_tran in matrix_ce_tran_file:
			# node object IDs to which the particular transmission line object is connected
			# node IDs of the node objects to which this transmission line is connected.
				ser_num = matrix_ce_tran['global_serial'] #global serial number of the shared existing transmission line
				t_node_id1 = matrix_ce_tran['from_node'] #From end node 
				node_zone1 = matrix_ce_tran['from_zone'] #From end zone number
				t_node_id2 = matrix_ce_tran['to_node'] #To end node
				node_zone2 = matrix_ce_tran['to_zone'] #To end zone number
			#Parameters for Transmission Line
			#reactance
				react = matrix_ce_tran['reactance']
			#values of maximum allowable power flow on line in the forward and reverse direction:
			#Forward direction:
				pt_max = matrix_ce_tran['line_limit']/100
				life_time = matrix_ce_tran['life_time'] #life time of the candidate line
				interest_rate = matrix_ce_tran['interest_rate'] #interest rate of the investment 
				cost_per_cap = matrix_ce_tran['cost_per_cap']*pt_max #capital cost for the construction 
				pres_absence = matrix_ce_tran['pres_absence'] #status of the construction 
				if ser_num in se_line_ser_list:
					new_node_base = 0
					for aggr_count in node_zone1:
						new_node_base = new_node_base + node_num_vector[aggr_count]
					t_node_id1 = t_node_id1 + new_node_base
					new_node_base = 0
					for aggr_count in node_zone2:
						new_node_base = new_node_base + node_num_vector[aggr_count]
					t_node_id2 = t_node_id2 + new_node_base
					cand_line_instance = self.cand_line( k + 1, ser_num, self.node_object[ t_node_id1 - 1 ], self.node_object[ t_node_id2 - 1 ], pt_max, react, interest_rate, life_time, cost_per_cap, pres_absence )
					cand_line_object.append( cand_line_instance )		# pushes the trans_line_instance object into the vector
					cand_line_ser_list.append(ser_num)  
					univ_cand_line_num = univ_cand_line_num + 1
					k +=1
			matrix_ce_tran_file.close() #Close 
##End of my translation
	assign_prob()
	#end constructor
	def assign_prob(self):
		for f in range(self.count_of_scenarios):
			self.probability.append(1/self.count_of_scenarios)

	def milp_avg_hr_gurobi(self): #Function MILPAvgHRGUROBI() implements the Mixed Integer Linear Programming Unit Commitment Solver routine by calling GUROBI routines for average heat rate objective for Horizontal Coordination Investment decision making
		milp_avg_hr_central(matrix_net_file, matrix_gen_file, matrix_tran)
 		#Function MILP() ends
