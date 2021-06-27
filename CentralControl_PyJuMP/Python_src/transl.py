#Member functions for class transmissionLine
#include transmissionLine class definition from transl.h, _node class definition from node.h
from Python_src.node import _node
from Python_src.log import log

class transmissionLine(object):
	def __init__(self, id_of_transl, node_connt1, node_connt2, powert_max, reactance): #constructor begins
		self.transl_id = id_of_transl
		self.conn_nodet1_ptr = node_connt1
		self.conn_nodet2_ptr = node_connt2
		self.pt_max = powert_max
		self.react = reactance
		self.device_nature = 0
		#log.info("\nInitializing the parameters of the transmission line with _id: {}".format(self.transl_id))
		self.from_node = self.conn_nodet1_ptr.get_node_id()
		self.to_node = self.conn_nodet2_ptr.get_node_id()
		self.conn_nodet1_ptr.sett_conn(self.id_of_transl, 1, self.react, self.to_node) #increments the txr line connection variable to node 1
		self.conn_nodet2_ptr.sett_conn(self.id_of_transl, -1, self.react, self.from_node) #increments the txr line connection variable to node 2
		# constructor ends

	#def __del__(): #destructor
		#log.info("\nThe transmission line object having _id {} have been destroyed.\n".format(self.transl_id))
		#end of destructor

	def get_transl_id(self): #function get_transl_id begins
		return self.transl_id #returns the _id of the generator object
		#end of get_transl_id function

	def get_flow_limit(self): #function get_flow_limit begins
		return self.pt_max #returns the Maximum power flow limit
		# end of get_flow_limit function

	def get_transl_node_id1(self): #function get_gen_node_id begins
		return self.conn_nodet1_ptr.get_node_id() #returns the _id number of the node to which the generator object is connected 
		#end of get_gen_node_id function

	def get_transl_node_id2(self): #function get_gen_node_id begins
		return self.conn_nodet2_ptr.get_node_id() #returns the _id number of the node to which the generator object is connected 
		#end of get_gen_node_id function

	def get_reactance(self):
		return self.react
	

