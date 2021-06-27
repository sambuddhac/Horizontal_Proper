#Member functions for class Load.
#include definition from node
from Python_src.node import Node
from Python_src.log import log

class Load(object):
	def __init__(self, id_of_load, node_connl, scenario_count, load_p): #constructor begins
		self.load_id = id_of_load
		self.conn_nodel_ptr = node_connl
		self.number_of_scenarios = scenario_count
		self.device_nature = 0
		self.set_load_value(load_p) #Sets the load for each scenario
		self.conn_nodel_ptr.set_l_conn(id_of_load, self.pl) #increments the load connection variable to node
		#constructor ends

	#def __del__(): #destructor
		#log.info("\nThe load object having _id {} have been destroyed.\n".format(self.load_id ))
		#end of destructor

	def set_load_value(self, load_p): #Sets the load for each scenario
		for i in range(number_of_scenarios):
			self.pl.append(load_p[i])

	def get_load_id(self): #function get_load_id begins
		return self.load_id #returns the _id of the load object
		#end of get_load_id function

	def get_load_node_id(self): #function get_load_node_id begins
		return self.conn_nodel_ptr.get_node_id() #returns the _id number of the node to which the load object is connected
		#end of get_load_node_id function
