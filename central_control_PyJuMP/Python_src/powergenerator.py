#Include definition of Node class
from Python_src.node import Node

#constructor definition for the piecewise linear objective function
class Powergenerator(object):
	def __init__(self, _id, node_conng, inc_cost, no_load, max, min):
		self.gen_id = _id #Initializer list to initialize data members that don't need validity check
		self.conn_nodeg_ptr = node_conng
		self.no_load_cost = no_load
		self.incremental_cost = inc_cost #Initialize this data member to zero (unused for piecewise linear objective)
		self.set_gen_params_simple(max, min) #call the set function to perform validity check on parameter value ranges and assign the values 
		self.conn_nodeg_ptr.setg_conn(self.gen_id) #increments the generation connection variable to node
		#end of constructor


	def set_gen_params_simple(self, max, min): #set function to set the Powergenerator class data members min max limits
		self.p_max = max
		self.pmin = min 
		#end of set_gen_params function

	def __del__(): #destructor definition
		#log.info("\n_generator object {}  destroyed".format(self.gen_id))
		#end of destructor

	def get_gen_id(self): #returns the Powergenerator _id number
		return self.gen_id
		#end of get_gen_id function

	def get_p_max(self): #function get_p_max begins
		return self.p_max
		#get_p_max ends

	def get_pmin(self): #function get_pmin begins
		return self.pmin
		#get_p_max ends

	def get_line_coeff(self): #Gets the linear coefficient (Incremental production cost
		return self.incremental_cost
		#function getLinCoeff ends

	def get_nl_cost(self): #Gets the no load cost
		return self.no_load_cost
		#function getNLCost ends 

