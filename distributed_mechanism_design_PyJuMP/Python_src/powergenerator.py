import numpy as np # typically this is what you want, and then instantiate numpy arrays wherever needed as np.array(). However, not sure, in this particular class, if we will need one
#import powergenerator #I am not sure abot this line #You don't need this line, since in Python (unlike C++), we don't create a separate header file for the class declaration
# Include definition of Node class
from Python_src.node import node
from math import *

# constructor definition for the linear objective function
class Powergenerator(object):
	def __init__(self, _id, node_conng, inc_cost, no_load, max, min): #Since Python does dynamic type checking, we don't typically need to specify the type of the variables
 		self.gen_id=_id # Initializer list to initialize data members that don't need validity check
 		self.conn_nodeg_ptr= node_conng
 		self.no_load_cost=no_load
 		self.incremental_cost = inc_cost # Initialize this data member to zero (unused for piecewise linear objective)
		self.set_gen_params_simple(max, min) # call the set_ function to perform validity check on parameter value ranges and assign the values
		self.conn_nodeg_ptr.set_g_conn(self.gen_id) # increments the generation connection variable to node
 	# end of constructor

	####### Make sure you take care of proper indentation. I have taken care of indentation here, if functions belong in the same module or class, you need to indent all those properlyarly
 	def set_gen_params_simple(self, max, min): # set_ function to set_ the Powergenerator class data members min max limits
		self.pmax = max
		self.pmin = min
 	# end of set_GenParams function
	
	def get_gen_id(self): # returns the Powergenerator _id number
		return self.gen_id
 	# end of get_gen_id function

	def get_pmax(self): # function get_pmax begins ####### change this function name to "get_pmax" instead of "Powergget_pmax"
		return self.pmax 
	# get_pmax ends

	def get_pmin(self): # function get_pmin begins
		return self.pmin
	# get_pmax ends

	def get_lin_Coeff(self): # Gets the linear coefficient (Incremental production cost
		return self.incremental_cost
	# function get_lin_Coeff ends

	def get_nl_cost(self): # Gets the no load cost 
		return self.no_load_cost
	# function get_nl_cost ends 
