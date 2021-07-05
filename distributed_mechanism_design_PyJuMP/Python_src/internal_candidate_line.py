#Member functions for class intCandLine.
from math import *
import numpy as np
#include <iostream>
#include <iomanip>
#include <cmath>
# include candLine class definition from intcandidateLine.h, Node class definition from node.h
#include "intcandidateLine.h"
from Python_src.node import node
#include "node.h"
#using namespace std;

class intCandLine(object):
	def _init_(self, id_of_transl, node_connt1, node_connt2, powert_max, reactance, roi, life, cap, abs_pres): # constructor begins
		self.transl_id=id_of_transl
	  	self.conn_nodet_ptr1=node_connt1
	  	self.conn_nodet_ptr2=node_connt2
	  	self.pt_max=powert_max
	  	self.react=reactance
	  	self.res_from_stage1=abs_pres
	  	self.status_of_construction=abs_pres
		from_node= self.conn_nodet_ptr1.get_node_id() 
		to_node=self.conn_nodet_ptr2.get_node_id()    
		self.conn_nodet_ptr1.set_int_cand_conn(id_of_transl, 1, self.react, to_node, self.status_of_construction) # increments the txr line connection variable to node 1
		self.conn_nodet_ptr2.set_int_cand_conn(id_of_transl, -1, self.react, from_node, self.status_of_construction) # increments the txr line connection variable to node 2
		
	###set_tran_data(self,cap, life, roi)
		# calls set_tran_data member function to set_ the parameter values

		# constructor ends

	def modify_node_react(self): # function to modify the nodal connected reactance, if the candidate line is actually built
		from_node=self.conn_nodet_ptr1.get_node_id()
		to_node=self.conn_nodet_ptr2.get_node_id()
		self.conn_nodet_ptr1.modify_react_app(self.transl_id, 1, self.react, to_node, 0) #increments the txr line connection variable to node 1
		self.conn_nodet_ptr2.modify_react_app(self.transl_id, -1, self.react, from_node, 0)  #increments the txr line connection variable to node 1
#	function ends

	def get_transl_id(self): # function get_transl_id begins
		return self.transl_id # returns the _id of the generator object
 # end of get_transl_id function


	def get_transl_node_id1(self): # function get_gen_node_id begins
		return self.conn_nodet_ptr1.get_node_id() # returns the _id number of the node to which the generator object is connected
# end of get_gen_node_id function

	def get_transl_node_id2(self): # function get_gen_node_id begins
		return self.conn_nodet_ptr2.get_node_id() # returns the _id number of the node to which the generator object is connected
# end of get_gen_node_id function

	def get_flow_limit(self): # Function get_flow_limit gets the value of power flow line limit	
		return self.pt_max
# Function get_flow_limit ends

	def get_reactance(self):
		return self.react	### I am not sure about self, actually I am a bit confused with "self"

	def set_tran_data(self, cap_cost, life_time, inter_rate): # member function to set_ parameter values of transmission lines
		self.capital_cost=cap_cost
		self.life_years=life_time
		self.rate_interest=inter_rate
# end function for set_ting parameter values

	def get_invest_cost(self): #member function get_invest_cost begins
		return (self.capital_cost*self.rate_interest*(pow((1+self.rate_interest), self.life_years)))/(pow((1+self.rate_interest), self.life_years)-1) #(1+rate_interest);capital_cost/100;//
	#return capital_cost;

	def return_pres_abs_status(self): # Returns the construction status of the candidate line
		return self.res_from_stage1

	def set_pres_abs_status(self): # Sets the construction status of the candidate line
		self.res_from_stage1=1 
