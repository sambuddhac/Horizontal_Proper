#Member functions for class Node.
from Python_src.log import log

class Node(object):
	def __init__(self, universal_id, id_of_node, zone_index): #constructor begins
		self.node_id = id_of_node #Node object id number
		self.zone_id = zone_index #Zone to which the node belongs
		#log.info("\nInitializing the parameters of the node with ID: {}".format(node_id))
		#initialize the connected devices to zero for node
		self.g_conn_number = 0 #number of generators connected to a particular node
		self.t_conn_number = 0 #number of transmission lines connected to a particular node
		self.l_conn_number = 0 #number of loads connected to a particular node
		self.shared_ex_conn_number = 0 #number of shared existing transmission lines connected to this node
		self.built_cand_conn_number = 0 # number of constructed candidate line connected to this node
		self.cand_conn_number = 0 # number of shared candidate transmission lines connected to this node
		self.int_cand_conn_number = 0 #number of internal candidate transmission lines connected to this node 
		self.shared_flag = 0 #node flag to indicate whether a shared existing or candidate line has been connected to a node #node flag to indicate whether a shared existing or candidate line has been connected to a node
		self.p_dev_count = 0 #initialize number of devices connectedto a node to zero #Number of devices connected to the particular node object
		self.from_react = 0.0 #Initialize the from reactance #Sum of reciprocals of reactances of lines for which this is the from node
		self.to_react = 0.0 #Initialize the to reactance #Sum of reciprocals of reactances of lines for which this is the to node
		self.global_rank = universal_id #sets the global_rank to universal_id #Global rank in the stitched list of shared line end nodes
		log.info("Node number: {}".format(global_rank))
		log.info("Zone number: {}".format(zone_id)) 
		self.gen_serial_num = [] #vector consisting of the serial numbers of generators connected to a particular node
		self.tran_from_serial = [] #vector consisting of the transmission lines for which the node is from node
		self.conn_node_list = [] #List of intra-zonal nodes that are directly connected to this node via transmission lines
		self.conn_react_rec = [] #List of reciprocals of reactances of the intra zone lines connected to the node
		self.tran_to_serial = []	#vector consisting of the transmission lines for which the node is to node
		self.conn_load_val = [] #Scenario values of connected load to this node
		self.load_serial_num = [] #vector consisting of the serial numbers of loads connected to a particular node
		self.se_from_serial = [] #vector consisting of the SE lines for which the node is from node
		self.se_to_serial = [] #vector consisting of the SE lines for which the node is to node
		self.share_node_list = [] #List of outer-zone nodes connected to the node via shared existing transmission lines
		self.share_react_rec = [] #List of reciprocals of reactances of the shared existing lines connected to the node
		self.cand_from_serial = [] #vector consisting of the cand lines for which the node is from node
		self.cand_to_serial = [] #vector consisting of the cand lines for which the node is to node
		self.int_cand_from_serial = [] #vector consisting of the internal cand lines for which the node is from node
		self.int_cand_to_serial = [] #vector consisting of the internal cand lines for which the node is to node
		self.built_cand_from_serial = [] #vector consisting of the built cand lines for which the node is from node
		self.built_cand_to_serial = [] #vector consisting of the built cand lines for which the node is to node	
		self.cand_node_list = [] #List of outer-zone nodes connected to the node via shared cand transmission lines
		self.cand_react_rec = [] #List of reciprocals of reactances of the shared cand lines connected to the node
		self.conn_shared_point = 0 #flag to indicate if this node is either the from or to end of any SE line. Default value 0 indicates it isn't
		self.conn_cand_point = 0 #flag to indicate if this node is either the from or to end of any Cand line. Default value 0 indicates it isn't
		self.conn_int_cand_point = 0 #flag to indicate if this node is either the from or to end of any internal Cand line. Default value 0 indicates it isn't 
		#constructor ends

	def __del__(self): #destructor
		#log.info("\nThe node object having ID {" <<  << "} have been destroyed.\n".format(self.node_id))
		#end of destructor

	def get_node_id(self): #function get_node_id begins
		return self.global_rank #returns node ID to the caller
		# end of function get_node_id

	def set_g_conn(self, serial_of_gen):
		self.g_conn_number += 1 #increment the number of generators connected by one whenever a generator is connected to the node
		self.gen_serial_num.append(serial_of_gen) #records the serial number of the generator connected to the node

	def set_t_conn(self, tran_id, dir, react, rank_of_other):
		self.t_conn_number += 1 #increment the number of txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.tran_from_serial.appendk(tran_id)
			self.from_react += 1/react	
			if rank_of_other in self.conn_node_list: #If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other) #find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] -= 1/react
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(-1/react)
		else:
			self.tran_to_serial.append(tran_id)
			self.to_react -= 1/react
			if rank_of_other in self.conn_node_list: #If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other) #find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] += 1/react
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(1/react)

	def setSEConn(self, tran_id, dir, react, rank_of_other):
		self.shared_ex_conn_number += 1 #increment the number of shared existing txr lines connected by one whenever a txr line is connected to the node
		self.conn_shared_point = 1 #set the flag to indicate this is the from or to end of a shared existing transmission line
		if dir == 1:
			self.se_from_serial.append(tran_id)
			self.from_react += 1/react
			if rank_of_other in self.conn_node_list: #If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other) #find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] -= 1/react
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(-1/react)
		else:
			self.se_to_serial.append(tran_id)
			self.to_react -= 1/react
			if rank_of_other in self.conn_node_list: #If predecided Gen value is given for this particular Powergenerator
				pos = conn_node_list.index(rank_of_other) #find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] += 1/react
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(1/react)
		log.info("The node {} in zone {} accessed by the SE line {} with reactance {} and dir {}".format(self.global_rank, self.zone_id, tran_id, react, dir))

	def setCandConn(self, tran_id, dir, react, rank_of_other):
		self.cand_conn_number += 1 #increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.cand_from_serial.append(tran_id)
		else:
			self.cand_to_serial.append(tran_id)
		self.conn_cand_point = 1 #Flag set to indicate that this node is connected to a cand line
		log.info("The node {} in zone {} accessed by the shared candidate line {} with reactance {}".format(self.node_id, self.zone_id, tran_id, react))

	def setIntCandConn(self, tran_id, dir, react, rank_of_other, const_stat):
		self.int_cand_conn_number += 1 #increment the number of shared cand txr lines connected by one whenever a txr line is connected to the node
		if dir == 1:
			self.int_cand_from_serial.append(tran_id)
			self.from_react += const_stat*(1/react)	
			if rank_of_other in self.conn_node_list: #If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other) #find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] -= const_stat*(1/react)
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(const_stat*(-1/react))
		else:
			self.int_cand_to_serial.append(tran_id)
			self.to_react -= const_stat*(1/react)
			if rank_of_other in self.conn_node_list: #If predecided Gen value is given for this particular Powergenerator
				pos = self.conn_node_list.index(rank_of_other) #find the position of the Powergenerator in the chart of predecided values
				self.conn_react_rec[pos] += const_stat*(1/react)
			else:
				self.conn_node_list.append(rank_of_other)
				self.conn_react_rec.append(const_stat*(1/react))
		self.conn_int_cand_point = 1 #Flag set to indicate that this node is connected to an internal cand line

	def set_l_conn(self, l_id, load_val):
		self.l_conn_number += 1 #increment the number of loads connected by one whenever a load is connected to the node
		self.load_serial_num.append(l_id)
		self.conn_load_val = []
		self.conn_load_val = load_val #total connected load

	def get_gen_length(self): #function get_node_id begins
		return len(self.gen_serial_num) #returns node ID to the caller 
		#end of function get_node_id

	def get_gen_ser(col_count):
		return self.gen_serial_num[col_count-1]

	#function redContNodeCount begins

	def init_load(self, scen_num): #Initialize the default loads on all nodes to zero
		for i in range(scen_num):
			self.conn_load_val.append(0)

	def devpinit_message(self, scen_c): #function devpinit_message begins
		return self.conn_load_val[scen_c] #return the total connected load 
		#function devpinit_message ends

	def get_shared_flag(self):
		return self.conn_shared_point #return the status if this node is connected to a shared existing line

	def get_cand_flag(self):
		return self.conn_cand_point #return the status if this node is connected to a shared cand line

	def get_to_react(self):
		return self.to_react #return the total reciprocal of reactances for which this is the to node

	def get_from_react(self):
		return self.from_react #return the total reciprocal of reactances for which this is the from node

	def get_con_node_length(self):
		return len(self.conn_node_list) #returns the length of the vector containing the connected intra-zonal nodes

	def get_conn_ser(self, col_count):
		return self.conn_node_list[col_count-1] #returns the serial number of the connected internal node at this position

	def get_conn_react(self, col_count):
		return self.conn_react_rec[col_count-1] #returns the serial number of the connected internal node at this position

	def get_cand_line_length_f(self):
		return len(self.cand_from_serial) #returns the number of cand lines connected to this from node

	def get_cand_line_length_t(self):
		return len(self.cand_to_serial) #returns the number of cand lines connected to this to node

	def get_cand_ser_f(self, col_count):
		return self.cand_from_serial[col_count-1] #returns the serial number of the cand line at this position

	def get_cand_ser_t(self, col_count):
		return self.cand_to_serial[col_count-1] #returns the serial number of the cand line at this position

	def get_int_cand_line_length_f(self):
		return len(self.int_cand_from_serial) #returns the number of cand lines connected to this from node

	def get_int_cand_line_length_t(self):
		return len(self.int_cand_to_serial) #returns the number of cand lines connected to this to node

	def get_int_cand_ser_f(self, col_count):
		return self.int_cand_from_serial[col_count-1] #returns the serial number of the cand line at this position

	def get_int_cand_ser_t(self, col_count):
		return self.int_cand_to_serial[col_count-1] #returns the serial number of the cand line at this position
