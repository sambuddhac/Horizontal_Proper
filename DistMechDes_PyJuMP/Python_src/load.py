#Member functions for class Load.
#include definition from node
from Python_src.node import Node
from Python_src.log import log

class Load(object):
	def __init__(self, idOfLoad, nodeConnl, scenarioCount, Load_P): #constructor begins
		self.loadID = idOfLoad
		self.connNodelPtr = nodeConnl
		self.numberOfScenarios = scenarioCount
		self.deviceNature = 0
		self.setLoadValue(Load_P) #Sets the load for each scenario
		self.connNodelPtr.setlConn(idOfLoad, self.Pl) #increments the load connection variable to node
		#constructor ends

	#def __del__(): #destructor
		#log.info("\nThe load object having ID {} have been destroyed.\n".format(self.loadID ))
		#end of destructor

	def setLoadValue(self, Load_P): #Sets the load for each scenario
		for i in range(numberOfScenarios):
			self.Pl.append(Load_P[i])

	def getLoadID(self): #function getLoadID begins
		return self.loadID #returns the ID of the load object
		#end of getLoadID function

	def getLoadNodeID(self): #function getLoadNodeID begins
		return self.connNodelPtr.getNodeID() #returns the ID number of the node to which the load object is connected
		#end of getLoadNodeID function
	
