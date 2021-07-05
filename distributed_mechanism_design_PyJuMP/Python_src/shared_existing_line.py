#Member functions for class SELine.
#include Node class definition from node.py
from Python_src.node import Node
from Python_src.log import log


class SELine(object):
	def __init__(self, shared_rank, id_of_transl, node_connt, from_n, from_z, to_n, to_z, zonal_id_num, powert_max, reactance): #constructor begins
		self.transl_id = id_of_transl
		self.shared_index = shared_rank
		self.conn_nodet_ptr = node_connt
		self.pt_max = powert_max
		self.react = reactance
		self.from_zone = from_z
		self.to_zone = to_z
		self.from_node = from_n
		self.to_node = to_n
		self.other_node_global = 0
		log.info("\nInitializing the parameters of the transmission line with _id: {}".format(self.transl_id))
		if self.from_zone==zonal_id_num:
			self.conn_nodet_ptr.set_se_conn(id_of_transl, 1, self.react, self.to_zone ) #increments the txr line connection variable to node 1
			self.from_to_flag=1
		else:
			self.conn_nodet_ptr.set_se_conn(id_of_transl, -1, self.react, self.from_zone ) #increments the txr line connection variable to node 1
			self.from_to_flag=-1
		#constructor ends

	def get_transl_id(self): #function get_transl_id begins
		return self.transl_id #returns the _id of the generator object 
		#end of get_transl_id function

	def get_intl_node_id(self): #function get_gen_node_id begins
		return self.conn_nodet_ptr.get_node_id() #returns the _id number of the node to which the generator object is connected
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

	def get_ext_node_rank(self): #function get_gen_node_id begins
		return self.other_node_rank #returns the _id number of the node to which the generator object is connected
		#end of get_gen_node_id function

	def get_ext_node_global_rank(self): #function get_ext_node_global_rank for the outside-zone node _id begins
		return self.other_node_global #returns the global rank of the outside-zone node to which the SE line object is connected
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
		self.conn_nodet_ptr.send_ext_node_info(self.other_node_rank, self.from_to_outer, self.react, 0)

	def get_reactance(self):
		return self.react

	def assign_rank(self, ranking): #assigns rank to the from/to node
		self.conn_nodet_ptr.assign_global_rank(ranking)

	def connect_rank(self, ranking): #assigns rank to other_node_global
		self.other_node_global = ranking
		self.conn_nodet_ptr.populate_global_conn(ranking)
