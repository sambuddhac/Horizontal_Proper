#Member functions for class intCandLine
#include class definition from intcandidateLine, _node class definition from node
from Python_src.node import _node
from Python_src.log import log

class intCandLine(object):
	def __init__(self, id_of_transl, node_connt1, node_connt2, powert_max, reactance, roi, life, cap, abs_pres ): #constructor begins
		self.transl_id = id_of_transl
		self.conn_nodet_ptr1 = node_connt1
		self.conn_nodet_ptr2 = node_connt2
		self.pt_max = powert_max
		self.react = reactance
		self.status = abs_pres
		self.from_node = self.conn_nodet_ptr1.get_node_id()
		self.to_node = self.conn_nodet_ptr2.get_node_id()
		self.conn_nodet_ptr1.set_int_cand_conn(id_of_transl, 1, self.react, self.to_node, self.status) #increments the txr line connection variable to node 1
		self.conn_nodet_ptr2.set_int_cand_conn(id_of_transl, -1, self.react, self.from_node, self.status) #increments the txr line connection variable to node 2
		self.set_tran_data(cap, life, roi) #calls set_tran_data member function to set the parameter values
		#constructor ends

	def __del__(): #destructor
		#log.info("\nThe transmission line object having _id {} have been destroyed.".format(self.transl_id))
		#end of destructor

	def get_transl_id(self): #function get_transl_id begins
		return self.transl_id #returns the _id of the generator object
		#end of get_transl_id function

	def get_transl_node_id1(self): #function get_gen_node_id begins
		return self.conn_nodet_ptr1.get_node_id() #returns the _id number of the node to which the generator object is connected
		#end of get_gen_node_id function

	def get_transl_node_id2(self): #function get_gen_node_id begins
		return self.conn_nodet_ptr2.get_node_id() #returns the _id number of the node to which the generator object is connected
		#end of get_gen_node_id function

	def get_flow_limit(self): #Function get_flow_limit gets the value of power flow line limit
		return self.pt_max
		#Function get_flow_limit ends

	def getreactance(self):
		return self.react

	def set_tran_data(self, cap_cost, life_time, inter_rate): #member function to set parameter values of transmission lines
		self.capital_cost = cap_cost
		self.life_years = life_time
		self.rate_interest = inter_rate
		#end function for setting parameter values

	def get_invest_cost(self): #member function get_invest_cost begins
		#return (self.capital_cost*self.rate_interest*(((1+self.rate_interest) ** self.life_years)))/(((1+self.rate_interest) ** self.life_years)-1) #1+self.rate_interest);self.capital_cost/100
		return self.capital_cost

	def return_pres_abs_status(self): #Returns the construction status of the candidate line
		return self.status

