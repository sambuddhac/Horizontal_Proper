#Member functions for class candLine
#include Node class definition
from Python_src.node import Node
from Python_src.log import log

class candLine(object):
	def __init__(self, shared_rank, id_of_transl, node_connt, from_n, from_z, to_n, to_z, zonal_id_num, powert_max, reactance, roi, life, cap, abs_pres, owner): #constructor begins
		self.transl_id = id_of_transl
		self.sharedIndex = shared_rank
		self.conn_nodet_ptr = node_connt
		self.pt_max = powert_max
		self.react = reactance
		self.from_zone = from_z
		self.to_zone = to_z
		self.from_node = from_n
		self.to_node = to_n
		self.other_node_global = 0
		self.res_from_stage1 = abs_pres
		self.global_rank = 0
		self.ownership = owner
		#log.info("\nInitializing the parameters of the transmission line with _id: {}".format(self.transl_id))
		if self.from_zone==zonal_id_num:
			self.conn_nodet_ptr.set__cand_conn(id_of_transl, 1, self.react, self.to_zone) #increments the txr line connection variable to node 1
			self.from_to_flag=1
		else:
			self.conn_nodet_ptr.set__cand_conn(id_of_transl, -1, self.react, self.from_zone) #increments the txr line connection variable to node 1
			self.from_to_flag=-1
		self.set_tran_data(cap, life, roi) #calls set_tran_data member function to set_ the parameter values

		#constructor ends

	def return_ownership(self): #Returns the value of ownership
		return self.ownership

	def modify_node_react(self): #function to modify the nodal connected reactance, if the candidate line is actually built
		#If the connected node is the from node
		if self.from_to_flag==1:
			self.conn_nodet_ptr.modify_react_app(self.transl_id, 1, self.react, self.other_node_rank, 1) #increments the txr line connection variable to node 1
			self.conn_nodet_ptr.send_ext_node_info(self.other_node_rank, self.from_to_outer, self.react, 0)
		#If the connected node is the to node
		else:
			self.conn_nodet_ptr.modify_react_app(self.transl_id, -1, self.react, self.other_node_rank, 1) #increments the txr line connection variable to node 1
			self.conn_nodet_ptr.send_ext_node_info(self.other_node_rank, self.from_to_outer, self.react, 0)
		#function ends

	def get_transl_id(self): #function get_transl_id begins
		return self.transl_id #returns the _id of the generator object
		#end of get_transl_id function

	def get_intl_node_id(self): #function get_gen_node_id for the intra-zone node _id begins
		return self.conn_nodet_ptr.get_node_id() #returns the _id number of the intra-zone node to which the candidate line object is connected
		#end of get_gen_node_id function

	def get_intl_zone_id(self): #returns _id number of intra-zonal node end zone to which the transmission line is connected
		if self.from_to_flag==1:
			return self.from_zone #returns the _id number of the from zone if the intra-zonal node is the from node
		else:
			return self.to_zone #returns the _id number of the to zone if the intra-zonal node is the to node 
	#end of get_intl_zone_id() function

	def get_ext_node_id(self): #returns _id number of outer-zonal node end to which the transmission line is connected
		if self.from_to_flag==1:
			return self.to_node #returns the _id number of the from zone if the intra-zonal node is the from node
		else: 
			return self.from_node #returns the _id number of the to zone if the intra-zonal node is the to node
	#end of get_gen_node_id function

	def get_ext_zone_id(self): #returns _id number of outer-zonal node end zone to which the transmission line is connected
		if self.from_to_flag==1:
			return self.to_zone #returns the _id number of the from zone if the intra-zonal node is the from node
		else:
			return self.from_zone #returns the _id number of the to zone if the intra-zonal node is the to node
	#end of get_gen_node_id function

	def get_ext_node_rank(self): #function get_gen_node_id for the outside-zone node _id begins
		return self.other_node_rank #returns the _id number of the outside-zone node to which the candidate line object is connected
	#end of get_gen_node_id function

	def get_ext_node_global_rank(self): #function get_ext_node_global_rank for the outside-zone node _id begins
		return self.other_node_global #returns the global rank of the outside-zone node to which the candidate line object is connected
	#end of get_ext_node_global_rank function

	def get_flow_limit(self): #Function get_flow_limit gets the value of power flow line limit
		return self.pt_max
	#Function get_flow_limit ends

	def get_flow_dir(self): #returns the value of the direction flag indicating whether the intra-zonal node end of the line is from (+1) or to (-1) end
		return self.from_to_flag
	#Function get_flow_dir ends

	def outer_node_index(self, rank_of_outer_node, dir_flag):
		self.other_node_rank=rank_of_outer_node
		self.from_to_outer=dir_flag
		self.conn_nodet_ptr.send_ext_node_info(self.other_node_rank, self.from_to_outer, self.react, 1)

	def get_reactance(self):
		return self.react

	def assign_rank(self, ranking): #assigns rank to the from/to node
		self.conn_nodet_ptr.assign_global_rank(ranking)

	def connect_rank(self, ranking): #assigns rank to other_node_global
		self.other_node_global = ranking
		self.conn_nodet_ptr.populate_global_conn(ranking)

	def get_other_zone(self): #function get_other_zone returns the _id number of the outside zone to which the other node is connected
		if self.from_to_flag==1:
			return self.to_zone #returns the _id number of the outside-zone to which the candidate line object is connected
		else:
			return self.from_zone
	#end of get_other_zone function

	def set_tran_data(self, cap_cost, life_time, inter_rate): #member function to set_ parameter values of transmission lines
		self.capital_cost=cap_cost
		self.life_years=life_time
		self.rate_interest=inter_rate
	#end function for set_ting parameter values

	def get_invest_cost(self): #member function get_invest_cost begins
		return (self.capital_cost*self.rate_interest*((1+self.rate_interest)**self.life_years))/(((1+self.rate_interest)**self.life_years)-1) #(1+rate_interest);capital_cost/100#
		#return capital_cost

	def assign_line_rank(self, glob_rank): #Assigns global rank to the candidate line
		self.global_rank = glob_rank #Global rank of the candidate line

	def get_global_rank(self): #Returns the global rank of the candidate line
		return self.global_rank #Global rank of the candidate line

	def return_pres_abs_status(self): #Returns the construction status of the candidate line
		return self.res_from_stage1

	def set_pres_abs_status(self): #Sets the construction status of the candidate line
		self.res_from_stage1=1
	#end
