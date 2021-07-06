#Member functions for class Node.
from math import *
import numpy as np
#include <iostream>
# include Node class definition from node.h
#include <vector>
#include <algorithm>
#include "node.h"
#using namespace std;

class Node(object):
	def _init_(self, id_of_node, zone_index): # constructor begins
		self.node_id=id_of_node
		self.zone_id=zone_index
#cout << "\nInitializing the parameters of the node with ID: " << node_id << endl;
#initialize the connected devices to zero for node
		self.g_conn_number = 0 # number of generators connected to a particular node
		self.t_conn_number = 0 # number of transmission lines connected to a particular node
		self.l_conn_number = 0 # number of loads connected to a particular node
		self.shared_ex_conn_number = 0  # number of shared existing transmission lines connected to this node
		self.built_cand_conn_number = 0 # number of constructed candidate line connected to this node
		self.built_int_cand_conn_number = 0 # number of constructed candidate line connected to this node
		self.cand_conn_number = 0 # number of shared candidate transmission lines connected to this node
		self.int_cand_conn_number = 0 # number of internal candidate transmission lines connected to this node 
		self.shared_flag = 0 # node flag to indicate whether a shared existing or candidate line has been connected to a node
		self.p_dev_count = 0 # initialize number of devices connectedto a node to zero
		self.from_react = 0.0 # Initialize the from reactance
		self.to_react = 0.0  # Initialize the to reactance
		self.global_rank = 0 # sets the global_rank to default value of 0 

	# constructor ends

	def get_node_id(self): # function get_node_id begins
		return self.node_id #returns node ID to the caller
 # end of function get_node_id

	def set_g_conn(self, serial_of_gen):
		self.g_conn_number+=1 # increment the number of generators connected by one whenever a generator is connected to the node
		self.gen_serial_num.append(serial_of_gen) # records the serial number of the generator connected to the node 
### 

	def set_t_conn(self, tran_id, dir, react, rank_of_other):
		self.t_conn_number+=1 # increment the number of txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.tran_from_serial.append(tran_id)
			self.from_react += 1/react	
			if  rank_of_other in self.conn_node_list:  # If predecided Gen value is given for this particular Powergenerator
				pos= self.conn_node_list.index(rank_of_other) # find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] -= 1/react
			else: 
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(-1/react)
	
		else:
			self.tran_to_serial.append(tran_id)
			self.to_react -= (1/react)
			if rank_of_other in self.conn_node_list: # If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other) # find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] += 1/react
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(1/react)

	def set_se_conn(self, tran_id, dir, react, connect_zone):
		self.shared_ex_conn_number +=1  # increment the number of shared existing txr lines connected by one whenever a txr line is connected to the node
		if  dir == 1:
			self.se_from_serial.append(tran_id)
			self.from_react += (1/react)		
		else:
			self.se_to_serial.append(tran_id)
			self.to_react -= (1/react)
		self.conn_shared_point=1 #Flag set to indicate that this node is connected to an SE line
		if connect_zone not in self.connected_zone_list: # If the connected zone isn't in the list
			self.connected_zone_list.append(connect_zone) # Put it on the list
			self.multiplicity +=1 # increase the multiplicity by 1

	def modify_react_app(self, tran_id, dir, react, rank_of_other, device_type): # Modifies the to and from reactances of lines connected to this node, to account for the newly constructed lines
		if device_type== 1:  # If shared candidate line
			self.built_cand_conn_number+=1 # increment the number of shared constructed candidate lines connected by one whenever the line is connected to the node
			if  dir == 1:  
				self.built_cand_from_serial.append(tran_id)
				self.from_react += (1/react)		
			else:
				self.built_cand_to_serial.append(tran_id)
				self.to_react -= (1/react)
	
		else: # If internal candidate line
			self.built_int_cand_conn_number+=1 # increment the number of shared constructed candidate lines connected by one whenever the line is connected to the node
			if  dir == 1:
				self.builtint_cand_from_serial.append(tran_id)
				self.from_react += (1/react)	
				if  rank_of_other in  self.conn_node_list: # If predecided Gen value is given for this particular Powergenerator
					pos = self.conn_node_list.index(rank_of_other)  # find the position of the Powergenerator in the chart of predecided values
					self.conn_react_rec[pos] -= (1/react)
				else:
					self.conn_node_list.append(rank_of_other)
					self.conn_react_rec.append((-1/react))
			else:
				self.builtint_cand_to_serial.append(tran_id)
				self.to_react -= (1/react)
				if  rank_of_other in self.conn_node_list: # If predecided Gen value is given for this particular Powergenerator
					pos = self.conn_node_list.index(rank_of_other) # find the position of the Powergenerator in the chart of predecided values
					self.conn_react_rec[pos] += (1/react)
				else:
					self.conn_node_list.append(rank_of_other)
					self.conn_react_rec.append((1/react))
		self.conn_built_cand_point=1

	def set_cand_conn(self, tran_id, dir, react, connect_zone):
		self.cand_conn_number+=1 # increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.cand_from_serial.append(tran_id)
		else:
			self.cand_to_serial.append(tran_id)
		self.conn_cand_point=1 # Flag set to indicate that this node is connected to a cand line
		if  connect_zone not in self.connected_zone_list: # If the connected zone isn't in the list
			self.connected_zone_list.append(connect_zone) # Put it on the list
			self.multiplicity+=1 # increase the multiplicity by 1

	def get_node_multiplicity(self): # get the multiplicity of the node i.e: the number of different zones (other than the one where it belongs) to which it is connected
		return self.multiplicity

	def set_int_cand_conn(self, tran_id, dir, react, rank_of_other, const_stat):
		self.int_cand_conn_number+=1 # increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.int_cand_from_serial.append(tran_id)
			self.from_react += const_stat*(1/react)	
			if  rank_of_other in self.conn_node_list: # If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other) # find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] -= const_stat*(1/react)
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(const_stat*(-1/react))
		else:
			self.int_cand_to_serial.append(tran_id)
			self.to_react -= const_stat*(1/react)
			if rank_of_other in self.conn_node_list: # If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other)  # find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] += const_stat*(1/react)
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(const_stat*(1/react))
		self.conn_int_cand_point=1	 # Flag set to indicate that this node is connected to an internal cand line

	def set_l_conn(self, l_id, load_val):
		self.l_conn_number +=1 # increment the number of loads connected by one whenever a load is connected to the node
		self.load_serial_num.append(l_id)
		self.conn_load_val = []
		self.conn_load_val = load_val ####

	def get_gen_length(self): # function get_node_id begins
		return self.gen_serial_num.len() # returns node ID to the caller 

	def get_gen_ser(self, col_count):
		return self.gen_serial_num[col_count-1] ###not sure at is

# function red_cont_node_count begins

	def init_load(self,scen_num): # Initialize the default loads on all nodes to zero
		i = 0 ### for (int i = 0; i < scen_num; ++i)
		for i in range(scen_num):
			self.conn_load_val.append(0) ###Not sure

	def devpinit_message(self,scen_c): # function devpinit_message begins
		return self.conn_load_val[scen_c] # return the total connected load ###Not sure
# function devpinit_message ends

	def send_ext_node_info(self, rank_of_outer, direction,reactance, indicator_se_cand): # Function to populate the connected outer-node list
		if rank_of_outer in self.share_node_list: # If outer node rank is present in the list
			pos = self.share_node_list.index(rank_of_outer)  # find the position of the outer node
			if direction == 1:
				self.share_react_rec[pos] -= (1-indicator_se_cand)*(1/reactance)
			else:
				self.share_react_rec[pos] += (1-indicator_se_cand)*(1/reactance)
		else:
			if direction == 1:
				self.share_node_list.append(rank_of_outer)
				self.share_react_rec.append((1-indicator_se_cand)*(-1/reactance))
			else:
				self.share_node_list.append(rank_of_outer)
				self.share_react_rec.append((1-indicator_se_cand)*(1/reactance))

	def get_shared_flag(self):
		return self.conn_shared_point # return the status if this node is connected to a shared existing line		

	def get_cand_flag(self):
		return self.conn_cand_point  # return the status if this node is connected to a shared cand line

	def get_built_cand_flag(self):
		return self.conn_built_cand_point #return the status if this node is connected to a shared cand line that is built

	def get_to_react(self):
		return self.to_react  #return the total reciprocal of reactances for which this is the to node

	def get_from_react(self):
		return self.from_react # return the total reciprocal of reactances for which this is the from node

	def get_con_node_length(self):
		return self.conn_node_list.len() # returns the length of the vector containing the connected intra-zonal nodes

	def get_conn_ser(self,col_count):
		return self.conn_node_list[col_count-1] # returns the serial number of the connected internal node at this position

	def get_conn_react(self,col_count):
		return self.conn_react_rec[col_count-1] # returns the serial number of the connected internal node at this position

	def get_extra_node_length(self):
		return self.share_node_list.len() # returns the length of the vector containing the connected outer-zonal nodes

	def get_ext_conn_ser(self,col_count):
		return self.share_node_list[col_count-1] # returns the serial number of the connected external node at this position

	def get_ext_conn_react(self,col_count):
		return self.share_react_rec[col_count-1] # returns the serial number of the connected internal node at this position

	def get_cand_line_length_f(self):
		return self.cand_from_serial.len() # returns the number of cand lines connected to this from node

	def get_cand_line_length_t(self):
		return self.cand_to_serial.len() # returns the number of cand lines connected to this to node

	def get_cand_ser_f(self,col_count):
		return self.cand_from_serial[col_count-1] #returns the serial number of the cand line at this position

	def get_cand_ser_t(self,col_count):
		return self.cand_to_serial[col_count-1]  # returns the serial number of the cand line at this position

	def get_int_cand_line_length_f(self):
		return self.int_cand_from_serial.len() # returns the number of cand lines connected to this from node

	def get_int_cand_line_length_t(self):
		return self.int_cand_to_serial.len()	# returns the number of cand lines connected to this to node

	def get_int_cand_ser_f(self,col_count):
		return self.int_cand_from_serial[col_count-1] # returns the serial number of the cand line at this position

	def get_int_cand_ser_t(self,col_count):
		return self.int_cand_to_serial[col_count-1]	# returns the serial number of the cand line at this position

	def assign_global_rank(self,rank): # Assigns the global rank to the nodes that are ends of shared lines
		self.global_rank = rank # sets the rank 
	def populate_global_conn(self,rank): # Populates the ext_node_global_rank vector with the global ranks
		if rank in self.share_node_list:	#If outer node rank is not present in the list
			self.ext_node_global_rank.append(rank)
	def get_global_rank(self): # Returns the global rank of this node
		return self.global_rank # Global rank in the stitched list of shared line end nodes
