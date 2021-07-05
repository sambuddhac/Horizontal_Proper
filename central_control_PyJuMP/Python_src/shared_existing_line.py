#Member functions for class SELine.
#include _node class definition from node.py
from Python_src.node import _node
from Python_src.log import log

class SELine(object):
	def __init__(self, local_rank, id_of_transl, node_connt1, node_connt2, powert_max, reactance): #constructor begins
		self.transl_id = id_of_transl
		self.local_index = local_rank
		self.conn_nodet_ptr1 = node_connt1
		self.conn_nodet_ptr2 = node_connt2
		self.pt_max = powert_max
		self.react = reactance	
		self.from_node = self.conn_nodet_ptr1.get__node_id()
		self.to_node = self.conn_nodet_ptr2.get__node_id()
		log.info("\nInitializing the parameters of the shared transmission line with _id: {}".format(self.transl_id))
		log.info("from node: {} To node: {}".format(self.from_node, self.to_node))
		self.conn_nodet_ptr1.set_se_conn(self.transl_id, 1, self.react, self.to_node) #increments the txr line connection variable to node 1
		self.conn_nodet_ptr2.set_se_conn(self.transl_id, -1, self.react, self.from_node) #increments the txr line connection variable to node 1 
		#constructor ends

	#def __del__(): #destructor
		#log.info("\nThe transmission line object having _id {} have been destroyed.\n".format(self.transl_id))
		# end of destructor

	def get_Transl_id(self): #function get_transl_id begins
		return self.transl_id #returns the _id of the generator object 
		#end of get_transl_id function

	def get_from_node_id(self): #function get_from_node_id begins
		return self.conn_nodet_ptr1.get_node_id() #returns the _id number of the from node 
		#end of get_from_node_id function

	def get_to_node_id(self): #function get_To_node_id begins
		return self.conn_nodet_ptr2.get_node_id() #returns the _id number of the to node
		# end of get_To_node_id function


	def get_from_limit(self): #Function get_from_limit get_s the value of power flow line limit
		return self.pt_max
		#Function get_from_limit ends

	def get_reactance(self):
		return self.react
