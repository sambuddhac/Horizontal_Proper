#Member functions for class candLine
#include Node class definition
from Python_src.node import Node
from Python_src.log import log

class candLine(object):
	def __init__(self, local_rank, id_of_transl, node_connt1, node_connt2, powert_max, reactance, rate_of_interest, life, cap, abs_pres): #constructor begins
		self.transl_id = id_of_transl
		self.local_index = local_rank
		self.conn_nodet_ptr1 = node_connt1
		self.conn_nodet_ptr2 = node_connt2
		self.pt_max = powert_max
		self.react = reactance
		self.status = abs_pres
		self.from_node = conn_nodet_ptr1.get_node_id()
		self.to_node = conn_nodet_ptr2.get_node_id()
		#log.info("\nInitializing the parameters of the transmission line with ID: {}".format(self.transl_id))
		self.conn_nodet_ptr1.set_cand_conn(id_of_transl, 1, self.react, self.to_node ) #increments the txr line connection variable to node 1
		self.conn_nodet_ptr2.set_cand_conn(id_of_transl, -1, self.react, self.from_node ) #increments the txr line connection variable to node 1
		self.set_tran_data(cap_cost, life, roi) #calls set_tran_data member function to set the parameter values
		#constructor ends

	#def __del__(): #destructor
		#log.info("\nThe transmission line object having ID {} have been destroyed.\n".format(self.transl_id))
		#end of destructor

	def get_transl_id(self): #function gettransl_id begins
		return self.transl_id #returns the ID of the generator object
		#end of gettransl_id function

	def get_from_node_id(self): #function get_from_node_id for the from node ID begins
		return self.conn_nodet_ptr1.get_node_id() #returns the ID number of the from  node 
		#end of get_from_node_id function

	def get_to_node_id(self): #function get_to_node_id for the to node ID begins
		return self.conn_nodet_ptr2.get_node_id() #returns the ID number of the to node 
		#end of get_to_node_id function

	def get_flow_limit(self): #Function get_flow_limit gets the value of power flow line limit
		return self.ptMax 
		#Function get_flow_limit ends

	def get_reactance(self):
		return self.react

	def set_tran_data(self, cap_cost, life_time, inter_rate): #member function to set parameter values of transmission lines
		self.capital_cost = cap_cost
		self.life_years = life_time
		self.rate_interest = inter_rate
		#end function for setting parameter values

	def get_invest_cost(self): #member function get_invest_cost begins
		#return (self.capital_cost*self.rate_interest*(((1+self.rate_interest) ** self.llife_years)))/(((1+self.rate_interest) ** self.llife_years)-1) #(1+self.rate_interest) #self.capital_cost/100
		return self.capital_cost

	def return_pres_abs_status(self): #Returns the construction status of the candidate line
		return self.status
